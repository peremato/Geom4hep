using Geom4hep
using CUDA
using CUDAKernels
using KernelAbstractions
using PrettyTables 

using BenchmarkTools, Test, Printf
#using BenchmarkPlots, StatsPlots

using Distributed
using ArgParse

#t = addprocs(nprocs())

const N = 1024*10

function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table s begin
        "--compare","-c"
            help = "compare with CPU version"
            action = :store_true
        "--cpugpu"
            help = "run CPU and GPU combined"
            action = :store_true           
        "--xray","-x"
            help = "run Xray benchmark"
            action = :store_true
        "--nprocs","-n"
            help = "number of processes"
            arg_type = Int
        "--nthreads","-t"
            help = "number of processes"
            arg_type = Int
         "--options","-o"
            help = "comma separated list of options cpu,gpu,mt,mp"
            arg_type = Any
    end

    return parse_args(s)
end


@everywhere begin

   using Geom4hep, Printf

   #---CPU implementation----------------------------------------------
   function cpu_d2out!(dist, shape, points, dirs, device = nothing)
	   n = length(points)
	   for i in 1:n
		   @inbounds dist[i] = distanceToOut(shape, points[i], dirs[i])
	   end
	   return dist
   end

   #---CPU parallel implementation----------------------------------------------
   function cpu_d2out_p!(dist, shape, points, dirs, device = nothing, from=1, to=length(points))
	   Threads.@threads for n in from:to
		 @inbounds dist[n] = distanceToOut(shape, points[n], dirs[n])
	   end
	   return dist[from:to]
   end

   function cpu_d2out_mp!(dist, shape, points, dirs, device = nothing)
	   nproc = nprocs();
	   block = ceil(Int, length(points)/nproc)
	   if nproc > 1   # Ensure at least one new worker is available
		   result = pmap(1:nproc) do n
				   from = (n-1)*block+1
				   to = n*block
				   return cpu_d2out_p!(dist, shape, points, dirs, device, from, to)
		   end
		   return(vcat(result...))
	   else
		 return (cpu_d2out_p!(dist, shape, points, dirs, device))
	   end
   end
end


 #---GPU implementation----------------------------------------------
 function cuda_kernel!(dist, shape, points, dirs)
	 n = length(points)
	 index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
	 stride = blockDim().x * gridDim().x
	 for i = index:stride:n
		 @inbounds dist[i] = distanceToOut(shape, points[i], dirs[i])
	 end
	 return
 end

 function cuda_d2out!(dist, shape, points, dirs, device = nothing)
	 dist_d = CUDA.fill(zero(typeof(dist[1])), N)
	 points_d = CuArray(points)
	 dirs_d = CuArray(dirs)
	 numblocks = ceil(Int, length(points)/256)
	 CUDA.@sync begin
		 @cuda threads=256 blocks=numblocks cuda_kernel!(dist_d, shape, points_d, dirs_d)
	 end
	 dist = Array(dist_d)
	 return dist
 end

 @kernel function d2out_kernel!(dist, shape, points, dirs)
	 N = length(points)
	 i = @index(Global)
	 @inbounds dist[i] = distanceToOut(shape, points[i], dirs[i])
 end

 function d2out!(dist,shape,points,dirs, device = "CPU")
	 n = device == "CUDA" ? 256 : 32
	 d = device == "CUDA" ? CUDADevice() : CPU();
	 dist_d   = device == "CUDA" ? CUDA.fill(zero(typeof(dist[1])), N) : dist
	 points_d = device == "CUDA" ? CuArray(points) : points
	 dirs_d   = device == "CUDA" ? CuArray(dirs) : dirs
	 numblocks = ceil(Int, length(points)/32)
	 kernel! = d2out_kernel!(d, n)
	 wait(kernel!(dist_d,shape,points_d,dirs_d,ndrange=size(points))) 
	 dist = Array(dist_d)
	 return (dist)
 end

function cpugpu!(n,dist,shape,points,dirs) 
 
   jobs    = Channel{Int}(n);
   results = Channel{Tuple}(n);

	function do_work(f,dist,shape,points,dirs,device)
	    job_id = take!(jobs)
	    result = f(dist,shape,points,dirs)
		put!(results, (job_id, device, 1))
	end;

	function make_jobs(n)
		 for i in 1:n
			 put!(jobs, i)
		 end
	end;

	errormonitor(@async make_jobs(n)); # feed the jobs channel with "n" jobs
	for i in 1:1
  	  errormonitor(@async do_work(cpu_d2out_p!,dist,shape,points,dirs,"CPU"))
    end
	for i in 1:1
  	   errormonitor(@async do_work(cuda_d2out!,dist,shape,points,dirs,"CUDA"))
    end     
    cpu = 0
    gpu = 0

    njobs = n
    
	@elapsed while n > 0 # print out results
		 job_id, device, result = take!(results)
		 if (device == "CPU")
		    cpu += 1
		    f = cpu_d2out_p! 
		 end
		 if (device == "CUDA")
		    gpu += 1
		    f = cuda_d2out!
		 end
		 n = n - 1
		 if n>0 
		   errormonitor(@async do_work(f,dist,shape,points,dirs,device))
		 end  
	end
	
	return (cpu,gpu)
	
end
 
function runBenchmarks(compare,cpugpu,options)

   d = Dict("cpu"  => (cpu_d2out!,"CPU"),
			"mt"   => (cpu_d2out_p!,"CPU"),
			"mp"   => (cpu_d2out_mp!,"CPU"),
			"gpu"  => (cuda_d2out!,"CUDA"),
			"kcpu" => (d2out!,"CPU"),
			"kgpu" => (d2out!,"CUDA")
   )

  results = []
    
  for T in (Float32, Float64)

    cone = Cone{T}(5, 10, 7, 15, 10, 0, 2π)
    tube = Tube{T}(5,10, 10, 0, 2π)
    box  = Box{T}(5, 10, 15)

    points = Vector{Point3{T}}(undef, N)
    dirs = Vector{Vector3{T}}(undef, N)
    
    for shape in (box, tube, cone)

        @printf "* Using shape: %s with %d points *\n" typeof(shape) N
        fillPoints!(points, shape, kInside)
        fillDirections!(dirs, shape, points)
        dist = fill(zero(T), N)

        if (compare)
           ref_dist = cpu_d2out!(dist, shape, points, dirs) 
        end
        
        for opt in options
            fcn,device = d[opt]
            @printf "%s (%s) : " uppercase(opt) "$device"  
	        dist = fcn(dist, shape, points, dirs, device)
            if (compare) 
		        @test ref_dist ≈ dist
		    end  
            t = @benchmark $fcn($dist, $shape, $points, $dirs, $device)
            b = @elapsed for i in 1:1000
              fcn(dist, shape, points, dirs, device)
            end
            display(t)
            @printf("\n")
		    @printf "** %s (%s) : (%f seconds)\n" uppercase(opt) "$device"  b 
            push!(results,[opt,device,typeof(shape),
						   minimum(t).time/1000,
						   maximum(t).time/1000,
						   median(t).time/1000,
						   mean(t).time/1000,
						   b,
						   mean(t).allocs,
						   mean(t).memory/1024])
        end
          
		if (cpugpu)
		   bm = @elapsed cpu,gpu = cpugpu!(1000, dist, shape, points, dirs) 
		   @printf "** CPU: %d + GPU %d : (%f seconds)\n" cpu gpu bm 
		end
		@printf "\n"
    end
  end
  printTable(results)
end  

function printTable(data)
   pretty_table(permutedims(hvcat(size(data,1), data...));
   header = (["Type","Device","Shape","Min","Max","Median","Mean","Elapsed","Allocs","Mem"]))
end

const idxs =  [2 3 1; 1 3 2; 1 2 3]
const rdxs =  [3 1 2; 1 3 2; 1 2 3]

function generateXRay(model::CuGeoModel, ::Type{T}, npoints::Number, view::Int=1; cuda::Bool=true) where T<:AbstractFloat
    world = model.volumes[1]
    lower, upper = extent(model.shapes[world.shapeIdx])
    idx = (idxs[view,1], idxs[view,2], idxs[view,3])
    rdx = (rdxs[view,1], rdxs[view,2], rdxs[view,3])
    ix, iy, iz = idx
    nx = ny = round(Int, sqrt(npoints))
    xrange = range(lower[ix], upper[ix], length = nx)
    yrange = range(lower[iy], upper[iy], length = ny)
    result = zeros(T, nx,ny)
    if cuda && CUDA.functional()
        threads = 256 # (8,8)
        blocks = cld.((nx,ny),threads)
        cu_result = CuArray(result)
        cu_model = cu(model)
        CUDA.@sync @cuda threads=threads blocks=blocks k_generateXRay(cu_result, cu_model, lower, upper, idx, rdx)
        return xrange, yrange, Array(cu_result)
    else
        generateXRay(result, model, lower, upper, idx, rdx)
        return xrange, yrange, result
    end
end


function  generateXRay(result::Matrix{T}, model::CuGeoModel, lower::Point3{T}, upper::Point3{T}, idx::Tuple, rdx::Tuple) where T<:AbstractFloat
    nx, ny = size(result)
    ix, iy, iz = idx
    xi, yi, zi = rdx
    px = (upper[ix]-lower[ix])/(nx-1)
    py = (upper[iy]-lower[iy])/(ny-1)
    _dir = (0, 0, 1)
    dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi]) 
    Threads.@threads  for i in 1:nx
        for j in 1:ny
        state = CuNavigatorState{T}(1)
        _point = (lower[ix]+ i*px, lower[iy]+ j*py, lower[iz]+ kTolerance(T))
        point = Point3{T}(_point[xi], _point[yi], _point[zi])
        locateGlobalPoint!( model, state, point)
        mass::T =  0.0
        step::T =  0.0
        while step >= 0.0
            vol = model.volumes[state.currentVol]
            density = model.materials[vol.materialIdx].density
            step = computeStep!(model, state, point, dir, T(1000.))
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
     end
   end 
end

function k_generateXRay(result, model, lower::Point3{T}, upper::Point3{T}, idx::Tuple, rdx::Tuple) where T<:AbstractFloat
    nx, ny = size(result)

    ix, iy, iz = idx
    xi, yi, zi = rdx

    px = (upper[ix]-lower[ix])/(nx-1)
    py = (upper[iy]-lower[iy])/(ny-1)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if i <= nx && j <= ny
        _dir = (0, 0, 1)
        dir = Vector3{T}(_dir[xi], _dir[yi], _dir[zi])

        _point = (lower[ix]+ i*px, lower[iy]+ j*py, lower[iz]+ kTolerance(T))
        point = Point3{T}(_point[xi], _point[yi], _point[zi])

        state = CuNavigatorState{T}(1)
        locateGlobalPoint!(model, state, point)
        mass::T =  0.0
        step::T =  0.0
        while step >= 0.0
            vol = model.volumes[state.currentVol]
            density = model.materials[vol.materialIdx].density
            step = computeStep!(model, state, point, dir, 1000.)
            point = point + dir * step
            mass += step * density
        end
        @inbounds result[i,j] = mass
    end
    return
end

function XRay(isCuda)
    T = Float32
    world = processGDML("examples/trackML.gdml",T)
    volume = world.daughters[2].volume.daughters[1].volume
    model = fillCuGeometry(getWorld(volume))
    #-----Generate maps----------------------------------------    
    @printf "generating x-projection\n"
    rx = generateXRay(model, T, 1e6, 1, cuda=isCuda)
    @printf "generating y-projection\n"
    ry = generateXRay(model, T, 1e6, 2, cuda=isCuda)
    @printf "generating z-projection\n"
    rz = generateXRay(model, T, 1e6, 3, cuda=isCuda)
end

function XRay()
	t2 = @elapsed XRay(false) 
	@show t2
	t1 = @elapsed XRay(true)
	@show t1
	@printf "speedup [CPU/GPU]= %f\n" t2/t1
end

function main()
   args = parse_commandline()
   compare = get!(args, "compare",  false)
   cpugpu = get!(args, "cpugpu",  false)
   xray = get!(args, "xray",  false)

   s = get!(args, "options","")
   if  s != nothing
	  options = split(s,",") 
	  for i in 1:length(options) 
		 options[i] = lowercase(options[i])
	  end
	  runBenchmarks(compare,cpugpu,options)
   end
   
   if (xray)
      XRay()
   end

end

main()
