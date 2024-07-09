module NyquistCircuitModeling

using Plots
using CSV
using NativeFileDialog
using DataFrames
using StatsBase
using DataInterpolations
using Optimization

export find_maxima_Zimg, Parameter_calculation, Estimate_Parameters, Theoretical_Z,
difference_et,optim_param

function plot_format(x,y,sc,xl,yl)
    max=maximum(x)
    min=minimum(x)
    plot(x,y,xscale=sc,
    framestyle=:box,right_margin=7*Plots.mm,linewidth=4,
    formatter=:plain,
    top_margin=5*Plots.mm,xlabel=xl,ylabel=yl,legend=false,
    #=ylims=(0,maximum(y)+maximum(y)/10)=#)
end

function Difference_plot(Xd,Yd,sc_d,xl_d,yl_d)
    D=diff(Yd)
    scatter(midpoints(Xd),D,xscale=sc_d,xlabel=xl_d,ylabel=yl_d,
    framestyle=:box,right_margin=7*Plots.mm,linewidth=4,
    legend=false,top_margin=5*Plots.mm)
end

function print_vector(v)
    if length(v) ==1
        (string(v))
    else
    join(string.(v),"; ")
    end
end

function find_maxima_Zimg(df)
    df = df[df."-Z'' (Ω)".>=0, :]
    Index=df."Index"
    Zimg=df."-Z'' (Ω)"
    M=[]
    I=[]
    m=[]
    ix=[]
    #print(maximum(Index))
    for i in 2:(length(Index)-1)
        #println(i," ",Zimg[i])
        if(Zimg[i-1]<Zimg[i]>Zimg[i+1])
            push!(M,Zimg[i])
            push!(I,i)
        else
            if (Zimg[i-1]>Zimg[i]<Zimg[i+1])
                push!(m,Zimg[i])
                push!(ix,i)
            end
        end
    end

    if (length(M)==1)
    println("Maxima = ",print_vector(M),"\n","Index = ",print_vector(I),"\n")
    printstyled("There is most likely $(length(M)) parallel RC circuit for this data file\n"
    ,bold=true)
    else
    println("Maxima = ",print_vector(M),"\n","Index = ",print_vector(I),"\n")
    printstyled("There are most likely $(length(M)) parallel RC circuits for this data file\n"
    ,bold=true)
    end

    if(length(m)==length(M))
        printstyled("Series RC's might be present in the circuit\n",bold=true)
        println("\n","minima = ",print_vector(m),"\n","indexm = ",print_vector(ix))
    else printstyled("Probably no series RC's are present in the circuit\n",bold=true)
        println("\n","minima = ",print_vector(m),"\n","indexm = ",print_vector(ix))
    end

    #print(m,ix)
    return M,I,m,ix
end

function Parameter_calculation(df,M,I,ix)
    f=df."Frequency (Hz)"
    Zre=df."Z' (Ω)"
    df = df[df."-Z'' (Ω)".>=0, :]
    R=[]
    C=[]
    Rp=[]
    Rs=Zre[1]
    push!(Rp,((Zre[I[1]])^2 + (M[1])^2)/Zre[I[1]]-Rs)
    #push!(Rp,R[1]-Rs)
    for i in 1:length(I)
        push!(R,((Zre[I[i]]).^2 .+ (M[i]).^2)./Zre[I[i]])
        push!(C,M[i]./((2*π*f[I[i]]).*(Zre[I[i]].^2 .+ M[i].^2)))
    end

    for i in 2:length(I)
        push!(Rp,R[i].-R[i-1])
    end
    println("\n","R=",R,"\n","C=",C,"\n","Rs=",Rs,"\n","Rp=",Rp)
end

function Parameter_calculation(df,M,I)
    f=df."Frequency (Hz)"
    Zre=df."Z' (Ω)"
    df = df[df."-Z'' (Ω)".>=0, :]
    R=[]
    C=[]
    Rp=[]
    Rs=Zre[1]
    push!(Rp,((Zre[I[1]])^2 + (M[1])^2)/Zre[I[1]]-Rs)
    for i in 1:length(I)
        push!(R,Zre[I[i]].*2)
    end

    for i in 2:length(I)
        push!(Rp,R[i].-R[i-1])
    end

    for i in 1:length(I)
        push!(C,1 ./(2*π*f[I[i]].*Rp[i]))
    end
    println("\n","R = ",print_vector(R),"\n","C = ",print_vector(C),"\n","Rs = ",
    print_vector(Rs),"\n","Rp = ",print_vector(Rp))
end

function Parameter_calculation(df,M,I,m,ix)
    f=df."Frequency (Hz)"
    Zre=df."Z' (Ω)"
    Index=df."Index"
    df = df[df."-Z'' (Ω)".>=0, :]
    R=[]
    C=[]
    Rp=[]
    Rs=Zre[1]
    
    for i in 1:length(ix)
        push!(R,Zre[ix[i]])
    end

    if(length(m)<length(M)) 
        push!(R,Zre[length(Index)])
    end

    push!(Rp,R[1]-Rs)

    if(length(I)>1)
        for i in 2:length(I)
        push!(Rp,R[i].-R[i-1])
        end
    end

    for i in 1:length(I)
        push!(C,1 ./(2*π*f[I[i]].*Rp[i]))
    end

    println("\n","R = ",print_vector(R),"\n","C = ",print_vector(C),"\n","Rs = ",print_vector(Rs),
    "\n","Rp = ",print_vector(Rp))

    return ix,Rs,Rp,C
end

function RC_series(df,ix,Zre)
    df = df[Zre .>= Zre[ix[length(ix)]], :]
    f=df."Frequency (Hz)"
    Zre=df."Z' (Ω)"
    Zimg=df."-Z'' (Ω)"
    Z=df."Z (Ω)"
    C=1 ./ (2*π*f .*Zimg)
    d=[]
    p_N_s = plot(midpoints(Zre),diff(C),legend=false)
end


function Estimate_Parameters()
    ff=pick_file()
    df=CSV.read(ff,DataFrame)
    df = df[df."-Z'' (Ω)".>=0, :] 

    Index=df."Index"
    f=df."Frequency (Hz)"
    Zre=df."Z' (Ω)"
    Zimg=df."-Z'' (Ω)"
    Z=df."Z (Ω)"
    Phase=df."-Phase (°)"
    Time="Time (s)"
    
    p_module=plot_format(f,Z,:log10,"Frequency (Hz)","Z (Ω)")
    savefig(p_module,ff*"_Module.html")

    p_N=plot_format(Zre,Zimg,:identity,"Zre (Ω)","Zimg (Ω)")
    savefig(p_N,ff*"_Nyquist.html")

    p_B=plot_format(f,Phase,:log10,"Frequency (Hz)","Phase Difference (deg)")
    savefig(p_B,ff*"_Bode.html")

    p_img=plot_format(f,Zimg,:log10,"Frequency (Hz)","Zimg (Ω)")
    savefig(p_img,ff*"_Zimg(f).html")
    
    p_ZimgDZre=plot_format(f,Zimg.-Zre,:log10,"Frequency (Hz)","Zimg-Zre (Ω)")
    savefig(p_ZimgDZre,ff*"_ZimgDZre.html")

    M,I,m,ix=find_maxima_Zimg(df)
    ix,Rs,Rp,C=Parameter_calculation(df,M,I,m,ix)

    RC_series(df,ix,Zre)
    return Rp,C,Rs
end


function Theoretical_Z(Rp,C,df)
    f=df."Frequency (Hz)"
    Zre_t= Rp[1] ./(1 .+ (2*π*f*Rp[1]*C[1]).^2) .+ Rp[2] ./(1 .+ (2*π*f*Rp[2]*C[2]).^2) .+
    Rp[3] ./(1 .+ (2*π*f*Rp[3]*C[3]).^2) .+ Rp[4] ./(1 .+ (2*π*f*Rp[4]*C[4]).^2)

    Zimg_t= (2*π*f.*(Rp[1].^2)*C[1]) ./ (1 .+ (2*π*f*Rp[1]*C[1]).^2) .+ 
    (2*π*f.*(Rp[2].^2)*C[2]) ./ (1 .+ (2*π*f*Rp[2]*C[2]).^2) .+
    (2*π*f.*(Rp[3].^2)*C[3]) ./ (1 .+ (2*π*f*Rp[3]*C[3]).^2) .+ 
    (2*π*f.*(Rp[4].^2)*C[4]) ./ (1 .+ (2*π*f*Rp[4]*C[4]).^2)

    return Zre_t,Zimg_t
end

function difference_et(x, df)
    Rp=x[1:length(x)÷2]
    C= x[length(x)÷2+1:length(x)]
    Zre_t,Zimg_t=Theoretical_Z(Rp,C,df)
    Zre=df."Z' (Ω)"
    Zimg=df."-Z'' (Ω)"
    sum(abs2.(Zre.-Zre_t))/length(Zre) + sum(abs2.(Zimg.-Zimg_t))/length(Zimg)
end

function optim_param()
    ff=pick_file()
    df=CSV.read(ff,DataFrame)
    M,I,m,ix=find_maxima_Zimg(df)
    ix,Rs,Rp,C=Parameter_calculation(df,M,I,m,ix)
    p0=Float64.(vcat(Rp,C))
    optimization_function=OptimizationFunction(difference_et,AutoForwardDiff())
    probleeem=OptimizationProblem(optimization_function,p0,df)
    rezz=solve(probleeem,Optimization.LBFGS(),maxiters=100)
    print(rezz-p0)
    return rezz
end

end
