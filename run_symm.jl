#default setting:
#filepath = "/home/xyz/xyz.w003.4000/xyz.w003/calculations/28712/"
filepath = "./"
const ev2Hartree = 0.036749308136649
const Bohr2Ang = 0.529177249

@info "initial_calculation config"
using JSON,DelimitedFiles,LinearAlgebra
str = open(io->read(io,String),string(filepath,"config"))
config = JSON.parse(str)
calc_job = config["calc_job"]

@info "environment setting"
if nworkers() == 1
    #push!(LOAD_PATH, dirname(dirname(dirname(@__FILE__))))
    Pkgpath = "/home/zounl/git/Hop_HT/Hop.jl"
    using Pkg
    Pkg.activate(Pkgpath)
    using Hop
else
    #@everywhere push!(LOAD_PATH, dirname(dirname(dirname(@__FILE__))))
    @everywhere Pkgpath = "/home/zounl/git/Hop_HT/Hop.jl"
    @everywhere using Pkg,JSON
    @everywhere Pkg.activate(Pkgpath)
    @everywhere using Hop
end


function genlist(x)
    return collect(range(x[1], stop = x[2], length = Int64(x[3])))
end

function k_data2kpath(kdata::AbstractString)
    return map(x->parse(Float64,x), split(kdata)[2:7])
end

function std_out_array(a::AbstractArray)
    return string(map(x->string(x," "),a)...)
end

@info "read TBmodel"
ts=time()
if config["interface"] == "openmx39"
    nm, fermi_level= Hop.Interface.createmodelopenmx(string(filepath,"openmx.scfout"))
elseif config["interface"] == "openmx38"
    nm, fermi_level= Hop.Interface.createmodelopenmx38(string(filepath,"openmx.scfout"))
elseif config["interface"] == "w90"
    nm=Hop.Interface.createmodelwannier(string(filepath,"w90.dat"))
    fermi_level = config["fermi_level"]
end
tf=time()
print("time for reading TBmodel is: ",tf-ts,"\r\n")

if config["is_symm"]
    ts=time()
    nm = convert(TBModel{ComplexF64}, nm)
    nm = changebasis(nm, Hop.Group.Us_openmx)
    set_is_canonical_ordered!(nm, true)
    
    @info "read symmtry operator & test if its a coset decomposition"
    symm =  Hop.Group.get_symmetries(string(filepath,"symm_O.json"),nm.lat,isdouble = nm.isspinful)
    group = Hop.Group.generategroup(symm)
    @assert length(symm) == length(group)
    if config["is_TR"]
        T = Hop.Group.getTRS(isdouble = nm.isspinful)
        push!(symm,T)
        group = Hop.Group.generategroup(symm)
    end

    @info "symmtrize & get_max_change"
    snm = Hop.Group.symmetrize(nm, group)
    print("max change of hamiltonian:\r\n")
    print(Hop.Group.get_max_change(snm,nm),"\r\n")
    nm = snm
    tf=time()
    print("time for symmtrization is: ",tf-ts,"\r\n")
end

if config["is_shared"]
    ts=time()
    @info "using shared_TB model"
    sm = SharedTBModel(nm)
    nm = sm
    tf=time()
    print("time for shared is: ",tf-ts,"\r\n")
end

if calc_job == "openmx_bs"  

    @info "get_kpath"
    str = open(io->read(io,String),string(filepath,"kpath.json"))
    k_data = JSON.parse(str)["k_path"]
    pnkpts = parse(Int64,split(k_data[1])[1])
    kpath = reshape(cat(dims=2,map(k_data2kpath,k_data)...),3,:)

    @info "calc_bandstructure"
    kdist, egvals = Hop.BandStructure.getbs(snm, kpath, pnkpts)

    @info "output_bandstructure in openmx format"
    f = open(string(filepath,"openmx.Band"),"w")
    println(f,snm.norbits," ",0," ",fermi_level)
    openmx_rlat = reshape((snm.rlat.*Bohr2Ang),1,:)
    println(f,std_out_array(openmx_rlat))
    println(f,length(k_data))
    for line in k_data
        println(f,line)
    end
    for i=1:div(size(kpath)[2],2)
        kstart = kpath[:,2*i-1]
        kend = kpath[:,2*i]
        k_list = zeros(Float64,pnkpts,3)
        for alpha = 1:3
            k_list[:,alpha] = genlist([kstart[alpha],kend[alpha],pnkpts])
        end
        for j = 1:pnkpts
            kvec = k_list[j,:]
            println(f,snm.norbits," ",std_out_array(kvec))
            println(f,std_out_array(ev2Hartree*egvals[:,(i-1)*pnkpts+j]))
        end
    end
    close(f)

elseif calc_job == "shiftcond_all"
    ωs = genlist(config["omegas"])
    @info "calculating shift conductivity"
    for alpha = 1:3, beta = 1:3, gamma = 1:3
        if config["isint"]
            sc = Hop.Optics.getshiftcond(nm,alpha,beta,gamma,ωs,
            fermi_level,convert(Array{Int64,1},
            config["kmesh"]),ϵ = config["epsilon"],isint=config["isint"], truncation=config["truncation"])
        else
            sc = Hop.Optics.getshiftcond(nm,alpha,beta,gamma,ωs,
            fermi_level,convert(Array{Int64,1},
            config["kmesh"]),ϵ = config["epsilon"])
        end
        open(string("shiftcond_",alpha,beta,gamma,".dat"), "w") do f
            writedlm(f, [ωs sc])
        end
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    @info "shiftcond_all calculated"
    
elseif calc_job == "shiftcond"
    ts=time()
    ωs = genlist(config["omegas"])
    @info "calculating shift conductivity"
    if config["isint"]
        @info "using intelligent kmesh"
        sc = Hop.Optics.getshiftcond(nm,config["alpha"],
        config["beta"],config["gamma"],ωs,
        fermi_level,convert(Array{Int64,1},
        config["kmesh"]),ϵ = config["epsilon"],isint=config["isint"], truncation=config["truncation"])
    else
        sc = Hop.Optics.getshiftcond(nm,config["alpha"],
        config["beta"],config["gamma"],ωs,fermi_level,
        convert(Array{Int64,1},
        config["kmesh"]),ϵ = config["epsilon"])
    end
    open(string("shiftcond_",config["alpha"],config["beta"],config["gamma"],".dat"), "w") do f
        writedlm(f, [ωs sc])
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    tf=time()
    print("time for shift conductivity calculation is: ",tf-ts,"\r\n")
    @info "shift conductivity calculated"

elseif calc_job == "SHG_all"
    ωs = genlist(config["omegas"])
    @info "calculating shift conductivity"
    for alpha = 1:3, beta = 1:3, gamma = 1:3
        σs = Hop.Optics.get_shg(nm,alpha,beta,gamma,ωs,
            fermi_level,convert(Array{Int64,1},
            config["kmesh"]),ϵ=config["epsilon"],scissor=config["scissor"])
        open(string("ReSHG_",alpha,beta,gamma,".dat"), "w") do f
            writedlm(f, [ωs real.(σs)])
        end
        open(string("ImSHG_",alpha,beta,gamma,".dat"), "w") do f
            writedlm(f, [ωs imag.(σs)])
        end
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    @info "SHG_all calculated"   

elseif calc_job == "SHG"
    ts=time()
    ωs = genlist(config["omegas"])
    @info "calculating shift conductivity"
    σs = Hop.Optics.get_shg(nm,config["alpha"],
        config["beta"],config["gamma"],ωs,fermi_level,
        convert(Array{Int64,1},
        config["kmesh"]),ϵ=config["epsilon"],scissor=config["scissor"])
    open(string("ReSHG_",config["alpha"],config["beta"],config["gamma"],".dat"), "w") do f
        writedlm(f, [ωs real.(σs)])
    end
    open(string("ImSHG_",config["alpha"],config["beta"],config["gamma"],".dat"), "w") do f
        writedlm(f, [ωs imag.(σs)])
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    tf=time()
    print("time for SHG calculation is: ",tf-ts,"\r\n")
    @info "SHG calculated"

elseif calc_job == "injection_all"
    ωs = genlist(config["omegas"])
    @info "calculating c current"
    for alpha = 1:3, beta = 1:3, gamma = 1:3
        cc = Hop.Optics.getinjection(nm,alpha,beta,gamma,ωs,
            fermi_level,convert(Array{Int64,1},
            config["kmesh"]),config["epsilon"])
        open(string("injection_",alpha,beta,gamma,".dat"), "w") do f
            writedlm(f, [ωs cc])
        end
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    @info "injection_all calculated"
    
elseif calc_job == "injection"
    ts=time()
    ωs = genlist(config["omegas"])
    @info "calculating injection current"
    cc = Hop.Optics.getinjection(nm,config["alpha"],
        config["beta"],config["gamma"],ωs,fermi_level,
        convert(Array{Int64,1},
        config["kmesh"]),config["epsilon"])
    open(string("injection_",config["alpha"],config["beta"],config["gamma"],".dat"), "w") do f
        writedlm(f, [ωs cc])
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    tf=time()
    print("time for injection calculation is: ",tf-ts,"\r\n")
    @info "injection calculated"

elseif calc_job == "jdos"
    ωs = genlist(config["omegas"])
    @info "calculating jdos"
    jdos = Hop.BandStructure.getjdos(nm, ωs, fermi_level,convert(Array{Int64,1},config["kmesh"]); ϵ=config["epsilon"])
    open(string("jdos.dat"), "w") do f
        writedlm(f, [ωs jdos])
    end
    open(string("result.dat"), "w") do f
        write(f, "calculation finish")
    end
    @info "jdos calculated"

elseif calc_job == "xyzHT"
    ts=time()
    @info "calculating xyzHT"
    @info "generate omegas"
    str = open(io->read(io,String),string(filepath,"openmx.Band.json"))
    bandjson = JSON.parse(str)
    gap = bandjson["band_gap"]
    @assert gap >= 0
    omega_low = min(gap-0.2,0)
    omega_up = 7
    omega_num = round(Int,(omega_up-omega_low)/config["omega_step"])
    omegas = [omega_low,omega_up,omega_num]
    print("HTomegas = ",omegas,"\r\n")
    ωs = genlist(omegas)
    @info "generate kpath"
    nkmesh = zeros(Int64,3)
    for i = 1:3
        nkmesh[i]=round(Int,config["kden"]*sqrt(dot(nm.rlat[:,i],nm.rlat[:,i]))/2)*2+1
    end
    @info "get tensor_dir"
    str = open(io->read(io,String),string(filepath,"tensor_dir.json"))
    tensor_dir = Dict{String,Array{Array{Int64,1},1}}(JSON.parse(str))
    σs = Hop.Xyz.get_xyz(nm, tensor_dir, ωs, fermi_level, nkmesh, config["epsilon"])
    idir=0
    for dir in tensor_dir["sc"]
        global idir
	idir = idir+1
        open(string("sc_",dir[1],dir[2],dir[3],".dat"), "w") do f
            writedlm(f, [ωs σs[:,idir]])
        end
    end
    for dir in tensor_dir["cc"]
        global idir
	idir = idir+1
        open(string("cc_",dir[1],dir[2],dir[3],".dat"), "w") do f
            writedlm(f, [ωs σs[:,idir]])
        end
    end
    for dir in tensor_dir["chi"]
        global idir
        idir = idir+1
        open(string("Re_chi_",dir[1],dir[2],".dat"), "w") do f
            writedlm(f, [ωs σs[:,idir]])
        end
        idir = idir+1
        open(string("Im_chi_",dir[1],dir[2],".dat"), "w") do f
            writedlm(f, [ωs σs[:,idir]])
        end
    end
    tf=time()
    print("time for xyz_HT calculation is: ",tf-ts,"\r\n")
    open(string("result.dat"), "w") do f 
        write(f, "calculation finish")
    end

end

