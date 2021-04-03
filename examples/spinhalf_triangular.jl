using LinearAlgebra
using QuantumHamiltonian
using LatticeTools
using Arpack

#=
    .   .   .   x   .   .   .

  .   .   x   o   o   x   .   .
               \ /
    .   .   .   O - o   .   .

  .   .   .   .   .   .   .   .

=#
function make_triangular_lattice(shape::AbstractMatrix{<:Integer})
    latticevectors = [1 -0.5; 0 0.5*sqrt(3.0)];
    unitcell = make_unitcell(latticevectors, SiteType=String)
    addsite!(unitcell, "A", carte2fract(unitcell, [0.0, 0.0]))
    nnbondtypes = [ [1, 0], [1, 1], [0, 1] ]
    nnnbondtypes = [ [ 2, 1], [ 1, 2], [-1, 1] ]

    lattice = make_lattice(unitcell, shape)
    hypercube = lattice.hypercube
    supercell = lattice.supercell
    tsym = FiniteTranslationSymmetry(lattice)
    psym = little_symmetry(tsym, PointSymmetryDatabase.find2d("6mm"))
    tsymbed = embed(lattice, tsym)
    psymbed = embed(lattice, psym)
    ssymbed = tsymbed ⋊ psymbed

    nnbonds = []
    nnnbonds = []

    for r_row in lattice.bravais_coordinates
        for colvec in nnbondtypes
            R_col, r_col = hypercube.wrap(r_row .+ colvec)
            roworb_super = ("A", r_row)
            colorb_super = ("A", r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(nnbonds, ((irow, icol), R_col))
        end
        for colvec in nnnbondtypes
            R_col, r_col = hypercube.wrap(r_row .+ colvec)
            roworb_super = ("A", r_row)
            colorb_super = ("A", r_col)
            irow = get(supercell.siteindices, roworb_super, -1)
            icol = get(supercell.siteindices, colorb_super, -1)
            push!(nnnbonds, ((irow, icol), R_col))
        end
    end

    return (unitcell=unitcell,
            lattice=lattice,
            space_symmetry_embedding=ssymbed,
            nearest_neighbor_bonds=nnbonds,
            next_nearest_neighbor_bonds=nnnbonds)
end

# shape = [  1  1; -1  2]  # 3
# shape = [  2  0;  0  2]  # 4
# shape = [  3  0;  0  3]  # 9
# shape = [  2  2; -2  4]  # 12
# shape = [  4  0;  0  4]  # 16
shape = [  3  3; -3  6]  # 27


trilat = make_triangular_lattice(shape)

n_sites = numsite(trilat.lattice.supercell)
@show n_sites

(hs, σ) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)
j1 = sum(σ(i1, j) * σ(i2 , j) for ((i1, i2), R) in trilat.nearest_neighbor_bonds for j in [:x, :y, :z])

H = 0.25 * j1

dimensions = Int[]
for (qn,) in quantum_number_sectors(hs)
    if qn < 0
        continue
    end
    @show qn
    hss = HilbertSpaceSector(hs, qn)
    hssr = represent_dict(hss)
    for ssic in get_irrep_components(trilat.space_symmetry_embedding)
        if qn == 0
            for z2_irrep in [
                [(GlobalBitFlip(false), 1), (GlobalBitFlip(true), 1)],
                [(GlobalBitFlip(false), 1), (GlobalBitFlip(true), -1)],
            ]
                irrep_iter = make_product_irrep(z2_irrep, get_irrep_iterator(ssic))
                rhssr = symmetry_reduce_parallel(hssr, irrep_iter)
                push!(dimensions, dimension(rhssr))
            end
        else
            irrep_iter = make_product_irrep(get_irrep_iterator(ssic))
            rhssr = symmetry_reduce_parallel(hssr, irrep_iter)
            push!(dimensions, dimension(rhssr))
            push!(dimensions, dimension(rhssr))
        end
    end
end

println("dimension:", dimensions)
println("total dimension: ", sum(dimensions))
println("max dimension: ", maximum(dimensions))
