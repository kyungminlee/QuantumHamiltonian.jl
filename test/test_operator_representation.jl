using Test
using QuantumHamiltonian

using Suppressor
using StaticArrays
using LinearAlgebra
using SparseArrays

using QuantumHamiltonian.Toolkit: pauli_matrix


@testset "OpRep" begin
    @testset "spin half" begin
        up = State("Up", +1)
        dn = State("Dn", -1)
        spinsite = Site([up, dn])
        n_sites = 4
        hs = HilbertSpace([spinsite for i in 1:n_sites])
        hsr = represent(hs)
        
        σ = Dict((isite, j) => pauli_matrix(hs, isite, j)
            for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);
            
        @testset "constructor" begin
            opr1 = OperatorRepresentation(hsr, σ[1, :x])
            opr2 = represent(hsr, σ[1, :x])
            @test opr1.hilbert_space_representation == hsr
            @test opr1.operator == σ[1, :x]
            @test opr2.hilbert_space_representation == hsr
            @test opr2.operator == σ[1, :x]
            @test opr1 == opr2
            @test get_space(opr1) === hsr
        end
            
        @testset "typetraits" begin
            opr1 = OperatorRepresentation(hsr, σ[1, :x])
            opr2 = OperatorRepresentation(hsr, σ[1, :y])
            OR = typeof(opr1)
            @test spacetype(OR) == typeof(hsr)
            @test operatortype(OR) == typeof(σ[1, :x])
            @test spacetype(opr1) == typeof(hsr)
            @test operatortype(opr1) == typeof(σ[1, :x])
            @test bintype(OR) == UInt
            @test bintype(opr1) == UInt
            
            @test scalartype(opr1) === Int
            @test scalartype(opr2) === Complex{Int}
            @test scalartype(typeof(opr1)) === Int
            @test scalartype(typeof(opr2)) === Complex{Int}
            
            @test valtype(opr1) === Int
            @test valtype(opr2) === Complex{Int}
            @test valtype(typeof(opr1)) === Int
            @test valtype(typeof(opr2)) === Complex{Int}
        end
            
        @testset "show" begin
            opr = OperatorRepresentation(hsr, σ[1, :x])
            show(devnull, MIME("text/plain"), opr)
        end
            
        @testset "properties" begin
            opr = OperatorRepresentation(hsr, σ[1, :x])
            dim = dimension(hsr)
            @test size(opr) == (dim, dim)
            @test size(opr, 1) == dim
            @test size(opr, 2) == dim
            @test_throws BoundsError size(opr, 3)
            
            @test bitwidth(opr) == n_sites
        end
            
        @testset "unary operator" begin
            op = σ[1, :x]
            opr = OperatorRepresentation(hsr, op)
            @test +opr == opr
            @test -opr == represent(hsr, -op)
        end
        
        @testset "binary operator" begin
            op1 = σ[1, :x]
            op2 = σ[3, :z]
            opr1 = OperatorRepresentation(hsr, op1)
            opr2 = OperatorRepresentation(hsr, op2)
            @test opr1 + opr2 == represent(hsr, op1 + op2)
            @test opr1 - opr2 == represent(hsr, op1 - op2)
            @test opr1 * opr2 == represent(hsr, op1 * op2)
            
            opr1 / 10
            @test opr1 * 2 == represent(hsr, op1 * 2)
            @test 2 * opr1 == represent(hsr, 2 * op1)
            
            opr1 / 10
            10 \ opr1
            op1 / 10
            10 \ op1
            
            @test opr1 / 2 == represent(hsr, op1 / 2)
            @test 2 \ opr1 == represent(hsr, 2 \ op1)
            @test opr1 // 2 == represent(hsr, op1 // 2)
            @test valtype(opr1//2) <: Rational{<:Integer}
            
            hs3 = HilbertSpace([spinsite, spinsite])
            op3 = pure_operator(hs3, 1, 1, 1)
            hsr3 = represent(hs3)
            opr3 = represent(hsr3, op3)
            @test_throws ArgumentError opr1 + opr3
            @test_throws ArgumentError opr1 - opr3
            @test_throws ArgumentError opr1 * opr3
        end
        
        @testset "symmetric" begin
            op1 = σ[2, :x]
            op2 = σ[3, :y]
            opr1 = OperatorRepresentation(hsr, op1)
            opr2 = OperatorRepresentation(hsr, op2)
            @test issymmetric(opr1)
            @test !issymmetric(opr2)
            @test ishermitian(opr1)
            @test ishermitian(opr2)
        end
        
        @testset "iterator" begin
            hs2 = HilbertSpace([spinsite, spinsite])
            op2 = pure_operator(hs2, 2, 1, 2)*2.0 + pure_operator(hs2, 2, 2, 1)*3.0 # 2 σ⁻₂ + 3 σ⁺₂
            hsr2 = represent(hs, UInt[0b00, 0b01])
            opr2 = represent(hsr2, op2)
            @test collect(get_row_iterator(opr2, 1)) == [-1 => 2.0]
            @test collect(get_row_iterator(opr2, 2)) == [-1 => 2.0]
            @test collect(get_column_iterator(opr2, 1)) == [-1 => 3.0]
            @test collect(get_column_iterator(opr2, 2)) == [-1 => 3.0]
        end
        
        @testset "get" begin
            opr = OperatorRepresentation(hsr, σ[2, :+]) # non-Hermitian
            dim = dimension(hsr)
            σ₊ = [0 1; 0 0]
            σ₀ = [1 0; 0 1]
            H0 = kron(σ₀, σ₀, σ₊, σ₀)
            opr_s = sparse(opr)
            opr_s1 = sparse_serial(opr)
            opr_s2 = sparse_parallel(opr)
            @test opr_s == opr_s1
            @test opr_s == opr_s2
            opr_d = Matrix(opr)
            @test opr_s == H0
            @test opr_d == H0
            @test opr[:,:] == opr_s
            @test isa(opr_s, SparseMatrixCSC)
            @test isa(opr_d, Matrix)
            for i in 1:dim
                @test get_row(opr, i) == H0[i, :]
                @test get_column(opr, i) == H0[:, i]
                for j in 1:dim
                    @test get_element(opr, i, j) == H0[i,j]
                    @test opr[i, j] == H0[i,j]
                end
                
                @test get_row(opr, i) != get_column(opr, i) # all rows and columns are different
                @test get_row(opr, i) == opr[i, :]
                @test get_column(opr, i) == opr[:, i]
            end
            @test_throws BoundsError get_row_iterator(opr, 0)
            @test_throws BoundsError get_row_iterator(opr, dim+1)
            @test_throws BoundsError get_column_iterator(opr, 0)
            @test_throws BoundsError get_column_iterator(opr, dim+1)
            @test_throws BoundsError get_row(opr, 0)
            @test_throws BoundsError get_row(opr, dim+1)
            @test_throws BoundsError get_column(opr, 0)
            @test_throws BoundsError get_column(opr, dim+1)
            @test_throws BoundsError get_element(opr, 0, 1)
            @test_throws BoundsError get_element(opr, dim+1, 1)
            @test_throws BoundsError get_element(opr, 1, 0)
            @test_throws BoundsError get_element(opr, 1, dim+1)
            @test_throws BoundsError opr[0, 1]
            @test_throws BoundsError opr[dim+1, 1]
            @test_throws BoundsError opr[1, 0]
            @test_throws BoundsError opr[1, dim+1]
        end # testset iterator
        
        @testset "apply" begin # complete space
            hsr_0 = represent(hs)
            dim = dimension(hsr_0)
            
            @testset "type" begin
                opr_x = represent(hsr_0, σ[2, :x])
                
                @test typeof(opr_x * rand(Int, dim)) === Vector{Int}
                @test typeof(opr_x * rand(Float64, dim)) === Vector{Float64}
                @test typeof(opr_x * rand(ComplexF64, dim)) === Vector{ComplexF64}
                
                opr_y = represent(hsr_0, σ[2, :y]) # operator has scalartype Complex{Int}
                @test typeof(opr_y * rand(Int, dim)) === Vector{Complex{Int}}
                @test typeof(opr_y * rand(Float64, dim)) === Vector{ComplexF64}
                @test typeof(opr_y * rand(ComplexF64, dim)) === Vector{ComplexF64}
            end
            
            @testset "apply!" begin
                for APP! in [apply!, apply_serial!, apply_parallel!]
                    opr = represent(hsr_0, σ[2, :+])
                    
                    σ₊ = [0 1; 0 0]
                    σ₀ = [1 0; 0 1]
                    
                    op_dense = kron(σ₀, σ₀, σ₊, σ₀)
                    state = rand(ComplexF64, dim)
                    
                    out1 = zeros(ComplexF64, dim)
                    
                    # matrix * columnvector
                    out0 = op_dense * state
                    APP!(out1, opr, state)
                    out2 = opr * state
                    @test isapprox(out0, out1, atol=1E-6)
                    @test isapprox(out0, out2, atol=1E-6)
                    
                    # add to the previous (do not overwrite)
                    APP!(out1, opr, state)
                    @test !isapprox(out0, out1, atol=1E-6)
                    
                    # rowvector * matrix
                    out1[:] .= zero(ComplexF64)
                    
                    out0 = transpose( transpose(state) * op_dense )
                    APP!(out1, state, opr)
                    out2 = state * opr
                    @test isapprox(out0, out1, atol=1E-6)
                    @test isapprox(out0, out2, atol=1E-6)
                    
                    # add to the previous (do not overwrite)
                    APP!(out1, state, opr)
                    @test !isapprox(out0, out1, atol=1E-6)
                    
                    let
                        state_long = rand(ComplexF64, dim+1)
                        @test_throws DimensionMismatch APP!(out1, state_long, opr)
                        @test_throws DimensionMismatch APP!(out1, opr, state_long)
                    end
                    
                    let
                        out_long = rand(ComplexF64, dim+1)
                        @test_throws DimensionMismatch APP!(out_long, state, opr)
                        @test_throws DimensionMismatch APP!(out_long, opr, state)
                    end
                end
            end
            
            @testset "mul!" begin
                opr = represent(hsr_0, σ[2, :x])
                vec_in = collect(1:dim) * 0.1
                
                let vec_out0, vec_out1
                    vec_out0 = zeros(Float64, dim)
                    vec_out1 = collect(1:dim) * 10.0
                    
                    apply!(vec_out0, opr, vec_in)
                    apply!(vec_out1, opr, vec_in)
                    
                    @test !isapprox(vec_out0, vec_out1, atol=1E-6)
                end
                
                let vec_out0, vec_out1
                    vec_out0 = zeros(Float64, dim)
                    vec_out1 = collect(1:dim) * 10.0
                    
                    mul!(vec_out0, opr, vec_in)
                    mul!(vec_out1, opr, vec_in)
                    
                    @test isapprox(vec_out0, vec_out1, atol=1E-6)
                end
                
            end
            # TODO(kyungminlee): Check for bounds error with range.
            
            @testset "error" begin
                hsr_small = represent(HilbertSpaceSector(hs, 0))
                opr = represent(hsr_small, σ[2, :x])
                dim_small = dimension(hsr_small)
                state = 2im*ones(ComplexF64, dim_small)
                
                out1 = zeros(ComplexF64, dim_small)
                tol = Base.rtoldefault(Float64)
                
                for APP! in [apply!, apply_serial!, apply_parallel!]
                    out1[:] .= zero(ComplexF64)
                    e1, e2 = APP!(out1, opr, state)
                    #@test isapprox(e1, 12im; atol=tol)
                    #@test isapprox(e2, 24.0; atol=tol)
                    @test all(isapprox.(out1, 0; atol=tol))
                    
                    e1, e2 = APP!(out1, state, opr)
                    #@test isapprox(e1, 12im; atol=tol)
                    #@test isapprox(e2, 24.0; atol=tol)
                    @test all(isapprox.(out1, 0; atol=tol))
                end
            end
            
        end
        
    end # testset spin half
end # testset OperatorRepresentation
