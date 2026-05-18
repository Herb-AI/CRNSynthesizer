using CRNSynthesizer
using Test, TestItems, TestItemRunner
using Aqua
using JET

#@testset "Aqua" Aqua.test_all(CRNSynthesizer)
#@testset "JET" JET.test_package(CRNSynthesizer; target_modules = [CRNSynthesizer])
@run_package_tests filter = ti -> ti.name == "Bottom-Up Molecule Synthesizer"
