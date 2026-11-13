using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Global (entire wing) sectional stiffness matrix of the Delft Pazy wing (nastran/Pazy_Delft/pazy_delft.bdf) from a static linear (SOL 101) analysis with NASTRAN. The wing is clamped at the root section and the six unit load cases are applied at the end of the reference axis.

reference_axis = 0.445

tips = [:flexible, :rigid]

for tip in tips
    T, L = nastran_pazy_delft_static_response(; reference_axis, tip)
    u_F1, u_F2, u_F3, u_M1, u_M2, u_M3 = (T[1:3, j] for j in 1:6)
    θ_F1, θ_F2, θ_F3, θ_M1, θ_M2, θ_M3 = (T[4:6, j] for j in 1:6)

    S = stiffness_matrix(L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

    write_stiffness_csv(pazy_delft_stiffness_file(; reference_axis, tip), S)

    # Spar alone (60 x 1.8 mm aluminum plate), as a reference
    b = 60e-3
    h = 1.8e-3
    E = 71e9
    ν = 0.33
    G = E / (2 * (1 + ν))
    A = b*h
    Iy = b*h^3/12
    Iz = h*b^3/12
    J = Iy+Iz
    Kt = 6.5e-3
    Ksy = Ksz = 5/6
    S0 = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz])

    println("Tip: $tip")
    println("Computed stiffness matrix:")
    display(round.(S, digits=4))
    println("Spar stiffness matrix:")
    display(S0)

    Smod = deepcopy(S)
    Smod[3,3] = 1e6

    uθ = displacements_from_loads(loadCase="all",L=L,S=S)
    uθMod = displacements_from_loads(loadCase="all",L=L,S=Smod)
    println("Generalized displacements from original stiffness matrix:")
    display(uθ)
    println("Generalized displacements from modified stiffness matrix:")
    display(uθMod)
    println("Difference in generalized displacements due to modified shear stiffnesses:")
    display(round.(uθMod .- uθ,digits=8))
    println("Ratio of modified to original generalized displacements:")
    display(round.(uθMod ./ uθ,digits=3))
end
