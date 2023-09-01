using PyCall 
using QuantumInformation
using LinearAlgebra
using Statistics
using Plots


symp = pyimport("Symplectic_generation")

# random unitary matrix 

gen=CUE(2)
# U=rand(gen)

n_qubit=2 # num qubit

total_run_number=10000

haar = HaarKet(2^n_qubit)
gen=CUE(2^n_qubit)

distance_clifford=zeros(total_run_number)
distance_unitary=zeros(total_run_number)
fidelity_clifford=zeros(total_run_number)
fidelity_unitary=zeros(total_run_number)
hs_clifford=zeros(total_run_number)
hs_unitary=zeros(total_run_number)

for l in 1:total_run_number

    ψ = rand(haar)

    
    ######### start random Clifford gate generation ##########
    Nmax_symp=convert(Int64, symp.numberofsymplectic(n_qubit)) # /4^n

    dim_sym=2^(n_qubit-1)
    randomSymplectic=symp.symplectic(rand(1:Nmax_symp),dim_sym)

    p=[sx,sy,sz,𝕀(2)]
    sel=rand(1:4,n_qubit)
    randomPauli=p[sel[1]]

    for i in 2:length(sel)   
        randomPauli=randomPauli⊗p[sel[i]] 
    end

    randomClifford=randomSymplectic*randomPauli

    ######### end random Clifford gate generation ##########


    # density matrix of state ψ
    σ = proj(ψ)   # initial

    # RUN Clifford
    
    #k = KrausOperators([sx, sy])
    k = KrausOperators([randomClifford])
    ρ = k(ψ)      # final
    ρ=ρ/tr(ρ)     #normalize

    #print(tr(ρ), "\t",tr(σ), "\n")
    
    distance_clifford[l]=trace_distance(ρ, σ)
    fidelity_clifford[l]=fidelity(ρ, σ)
    hs_clifford[l]=hs_distance(ρ, σ)
    
    # RUN Unitary  
    k = KrausOperators([rand(gen)])
    ρ = k(ψ)
    
    distance_unitary[l]=trace_distance(ρ, σ)
    fidelity_unitary[l]=fidelity(ρ, σ)
    hs_unitary[l]=hs_distance(ρ, σ)
end

plot(1:total_run_number, distance_clifford[1:total_run_number], label="clifford")
plot!(1:total_run_number, distance_unitary[1:total_run_number], label="unitary")


plot(1:total_run_number, fidelity_clifford[1:total_run_number], label="clifford")
plot!(1:total_run_number, fidelity_unitary[1:total_run_number], label="unitary")


plot(1:total_run_number, hs_clifford[1:total_run_number], label="clifford")
plot!(1:total_run_number, hs_unitary[1:total_run_number], label="unitary")
