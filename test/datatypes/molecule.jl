@testitem "BondType" begin
    @test instances(BondType) == (single, double, aromatic, triple, quadruple)

    @testset "to_string" begin
        @test to_string(single) == "-"
        @test to_string(double) == "="
        @test to_string(aromatic) == ":"
        @test to_string(triple) == "≡"
        @test to_string(quadruple) == "≣"
    end
end

@testitem "Atom" begin
    atom = Atom("H")
    @test atom.name == "H"

    atom = Atom("O")
    @test atom.name == "O"
end

@testitem "Bond" begin
    bond = Bond(1, 2, single)
    @test bond.from == 1
    @test bond.to == 2
    @test bond.bond_type == single

    bond = Bond(2, 1, double)
    @test bond.from == 1
    @test bond.to == 2
    @test bond.bond_type == double

    bond = Bond(3, 4, triple)
    @test bond.from == 3
    @test bond.to == 4
    @test bond.bond_type == triple
end

@testitem "Molecule" begin
    @testset "Constructor" begin
        atoms = [Atom("H"), Atom("O")]
        bonds = [Bond(1, 2, single)]
        molecule = Molecule(atoms, bonds, "[H]-[O]", Vector{UInt8}())

        @test length(molecule.atoms) == 2
        @test length(molecule.bonds) == 1
        @test molecule.atoms[1].name == "H"
        @test molecule.atoms[2].name == "O"
        @test molecule.bonds[1].from == 1
        @test molecule.bonds[1].to == 2
        @test molecule.bonds[1].bond_type == single
    end

    @testset "Count atoms" begin
        atoms = [Atom("H"), Atom("O"), Atom("H")]
        bonds = [Bond(1, 2, single), Bond(2, 3, single)]
        molecule = Molecule(atoms, bonds, "[H]-[O]-[H]", Vector{UInt8}())

        atom_counts = count_atoms(molecule)
        @test atom_counts["H"] == 2
        @test atom_counts["O"] == 1
    end

    @testset "From SMILES" begin
        smiles = "[C]-1(=[N]-[O]-1)-[H]"
        molecule = from_SMILES(smiles)
        @test length(molecule.atoms) == 4
        @test length(molecule.bonds) == 4
    end

    @testset "To SMILES" begin
        atoms = [Atom("[H]"), Atom("[O]"), Atom("[H]")]
        bonds = [Bond(1, 2, single), Bond(2, 3, single)]
        molecule = Molecule(atoms, bonds, "[H]-[O]-[H]", Vector{UInt8}())
        smiles = to_SMILES(molecule)
        @test smiles == "[H]-[O]-[H]"
    end

    @testset "SMILES roundtrip" begin
        smiles = "[H]-[H]"
        molecule = from_SMILES(smiles)
        smiles2 = to_SMILES(molecule)
        @test smiles == smiles2

        smiles = "[H]-[O]-[H]"
        molecule = from_SMILES(smiles)
        smiles2 = to_SMILES(molecule)
        @test smiles == smiles2

        smiles = "[O](-[O]-1)-[O]-1"
        molecule = from_SMILES(smiles)
        smiles2 = to_SMILES(molecule)
        @test smiles == smiles2

        smiles = "[C](-[H])(=[N]-1)-[O]-1"
        molecule = from_SMILES(smiles)
        smiles2 = to_SMILES(molecule)
        @test smiles == smiles2

        smiles = "[H]-[N](-[H])-[O]-[H]"
        molecule = from_SMILES(smiles)
        smiles2 = to_SMILES(molecule)
        @test smiles == smiles2
    end
end
