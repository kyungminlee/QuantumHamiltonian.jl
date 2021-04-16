using Test
using QuantumHamiltonian

using QuantumHamiltonian.Toolkit: pauli_matrix

@testset "Operator Application" begin
    
    @testset "spinhalf" begin
        QN = Tuple{Int}
        up = State("Up", 1)
        dn = State("Dn", -1)
        spin_site = Site([up, dn])
        hs = HilbertSpace([spin_site, spin_site, spin_site, spin_site])
        σ(i::Integer, j::Symbol) = pauli_matrix(hs, i, j)
        
        hsr = represent(HilbertSpaceSector(hs, 0))
        psi = SparseState{Float64, UInt}(UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
        
        @testset "apply" begin
            #psi = SparseState{Float64, UInt}(UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
            
            @test apply(psi, NullOperator()) == SparseState{Float64, UInt}()
            
            @test apply(psi, σ(1, :+)) == SparseState{Float64, UInt}()
            @test apply(psi, σ(1, :-)) == SparseState{Float64, UInt}(UInt(0b0010) => 2.0, UInt(0b0100) => 10.0)
            
            @test apply(psi, σ(2, :+)) == SparseState{Float64, UInt}(UInt(0b0111) => 10.0)
            @test apply(psi, σ(2, :-)) == SparseState{Float64, UInt}(UInt(0b0001) => 2.0)
            
            @test apply(NullOperator(), psi) == SparseState{Float64, UInt}()
            
            @test apply(σ(1, :+), psi) == SparseState{Float64, UInt}(UInt(0b0010) => 2.0, UInt(0b0100) => 10.0)
            @test apply(σ(1, :-), psi) == SparseState{Float64, UInt}()
            
            @test apply(σ(2, :+), psi) == SparseState{Float64, UInt}(UInt(0b0001) => 2.0)
            @test apply(σ(2, :-), psi) == SparseState{Float64, UInt}(UInt(0b0111) => 10.0)
            
            for i1 in 1:4, j1 in [:x, :y, :z, :+, :-]
                for i2 in 1:4, j2 in [:x, :y, :z, :+, :-]
                    let
                        ϕ1 = apply(apply(psi, σ(i1, j1)), σ(i2, j2))
                        ϕ2 = apply(psi, σ(i1, j1) * σ(i2,j2))
                        @test ϕ1 == ϕ2
                    end
                    let
                        ϕ1 = apply(σ(i1, j1), apply(σ(i2, j2), psi))
                        ϕ2 = apply(σ(i1, j1) * σ(i2,j2), psi)
                        @test ϕ1 == ϕ2
                    end
                end
            end
        end
        
        @testset "apply!" begin
            let
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, psi, NullOperator())
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, NullOperator(), psi)
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
            end
            
            let
                #psi = SparseState{Float64, UInt}(UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, psi, σ(1, :+))
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, σ(1, :+), psi)
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 2.25, UInt(0b0100) => 10.0)
            end
            
            let
                #psi = SparseState{Float64, UInt}(UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, psi, σ(1, :-))
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 2.25, UInt(0b0100) => 10.0)
                
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, σ(1, :-), psi)
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
            end
            
            let
                #psi = SparseState{Float64, UInt}(UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, psi, σ(2, :+))
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25, UInt(0b0111) => 10.0)
                
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, σ(2, :+), psi)
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25, UInt(0b0001) => 2.0)
            end
            
            let
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, psi, σ(2, :-))
                @test psi2 == SparseState{Float64, UInt}(UInt(0b0010) => 0.25, UInt(0b0001) => 2.0)
                
                psi2 = SparseState{Float64, UInt}(UInt(0b0010) => 0.25)
                apply!(psi2, σ(2, :-), psi)
                @test psi2 == SparseState{Float64, UInt}( UInt(0b0010) => 0.25, UInt(0b0111) => 10.0)
            end
            
            for i1 in 1:4, j1 in [:x, :z]
                ϕ1 = SparseState{ComplexF64, UInt}()
                ϕ2 = SparseState{ComplexF64, UInt}()
                apply!(ϕ1, psi, σ(i1, j1))
                apply!(ϕ2, σ(i1, j1), psi)
                @test ϕ1 == ϕ2
            end
            
            for i1 in 1:4, j1 in [:x, :y, :z, :+, :-]
                ϕ1 = apply(psi, σ(i1, j1))
                ϕ2 = SparseState{ComplexF64, UInt}()
                apply!(ϕ2, psi, σ(i1, j1))
                @test ϕ1 == ϕ2
            end
            
            let
                psi2 = SparseState{ComplexF64, UInt}(UInt(0b0010) => 0.25)
                psi3 = apply(psi2, σ(2, :y))
                psi4 = SparseState{Float64, UInt}()
                @test_throws InexactError apply!(psi4, psi2, σ(2, :y))
                # TODO more
            end
        end # testset apply!
        
        
        # @testset "materialize" begin
            #   n = 4
            #   heisenberg = sum(2 * σ(i, :+) * σ( mod(i, n) + 1, :-) for i in 1:n)
            #   heisenberg += sum(2 * σ(i, :-) * σ( mod(i, n) + 1, :+) for i in 1:n)
            #   heisenberg += sum(σ(i, :z) * σ( mod(i, n) + 1, :z) for i in 1:n)
            #   field_x = sum(σ(i, :x) for i in 1:n)
            #   field_y = sum(σ(i, :y) for i in 1:n)
            #   field_z = sum(σ(i, :z) for i in 1:n)
            #   @testset "zerosector" begin
                #     hsr = represent(HilbertSpaceSector(hs, 0))
                #     H, ϵ = materialize(hsr, heisenberg)
                #     @test H ≈ [ 0  2  0  0  2  0;
                #                 2 -4  2  2  0  2;
                #                 0  2  0  0  2  0;
                #                 0  2  0  0  2  0;
                #                 2  0  2  2 -4  2;
                #                 0  2  0  0  2  0]
                #     Hp, ϵp = materialize_parallel(hsr, heisenberg)
                #     @test Hp ≈ H
                #
                #     heisenberg2 = sum(σ(i, j) * σ(mod(i,n)+1, j) for i in 1:n, j in [:x, :y, :z])
                #     H2, ϵ2 = materialize(hsr, heisenberg2)
                #     @test H2 ≈ H
                #     H3, ϵ3 = materialize_parallel(hsr, heisenberg2)
                #     @test H3 ≈ H
                #
                #     let
                #       A, ϵ = materialize(hsr, NullOperator())
                #       Ap, ϵp = materialize_parallel(hsr, NullOperator())
                #       n = dimension(hsr)
                #       @test isempty(A.nzval)
                #       @test size(A) == (n, n)
                #       @test isempty(Ap.nzval)
                #       @test size(Ap) == (n, n)
                #     end
                #
                #     let
                #       Bx, ϵx = materialize(hsr, field_x)
                #       By, ϵy = materialize(hsr, field_y)
                #       Bz, ϵz = materialize(hsr, field_z)
                #       Bxp, ϵxp = materialize_parallel(hsr, field_x)
                #       Byp, ϵyp = materialize_parallel(hsr, field_y)
                #       Bzp, ϵzp = materialize_parallel(hsr, field_z)
                #       @test ϵx > 0
                #       @test ϵy > 0
                #       @test ϵz ≈ 0
                #       @test ϵx ≈ ϵxp
                #       @test ϵy ≈ ϵyp
                #       @test ϵz ≈ ϵzp
                #     end
                #   end
                #
                #   @testset "allsector" begin
                    #     PAULI_MATRIX = [
                    #       [0.0 1.0; 1.0 0.0], [0.0 -1.0im; 1.0im 0.0], [1.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 1.0]
                    #     ]
                    #     hsr = represent(hs)
                    #     H, ϵ = materialize(hsr, heisenberg)
                    #     H2  = sum( kron(PAULI_MATRIX[i], PAULI_MATRIX[i], PAULI_MATRIX[4], PAULI_MATRIX[4]) for i in 1:3)
                    #     H2 += sum( kron(PAULI_MATRIX[4], PAULI_MATRIX[i], PAULI_MATRIX[i], PAULI_MATRIX[4]) for i in 1:3)
                    #     H2 += sum( kron(PAULI_MATRIX[4], PAULI_MATRIX[4], PAULI_MATRIX[i], PAULI_MATRIX[i]) for i in 1:3)
                    #     H2 += sum( kron(PAULI_MATRIX[i], PAULI_MATRIX[4], PAULI_MATRIX[4], PAULI_MATRIX[i]) for i in 1:3)
                    #     @test H ≈ H2
                    #
                    #     let
                    #       Bx, _ = materialize(hsr, field_x)
                    #       By, _ = materialize(hsr, field_y)
                    #       Bz, _ = materialize(hsr, field_z)
                    #       @test isa(Bx.nzval, Array{Float64, 1})
                    #       @test isa(By.nzval, Array{ComplexF64, 1})
                    #       @test isa(Bz.nzval, Array{Float64, 1})
                    #     end
                    #     let
                    #       Bx, _ = materialize_parallel(hsr, field_x)
                    #       By, _ = materialize_parallel(hsr, field_y)
                    #       Bz, _ = materialize_parallel(hsr, field_z)
                    #       @test isa(Bx.nzval, Array{Float64, 1})
                    #       @test isa(By.nzval, Array{ComplexF64, 1})
                    #       @test isa(Bz.nzval, Array{Float64, 1})
                    #     end
                    #   end
                    # end
                end
            end
            