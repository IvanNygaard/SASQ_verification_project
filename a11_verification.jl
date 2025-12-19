# -----------------------------------------------------
# File: a11_verification.jl
# Author: Ivan Nygaard
# Created: 28/10/2025
# Description: Code used to verify the derivation of
#              the first element of the CC-Jacobian
#              that was derived by hand for the simplified
#              case of T2 = 0. 
# -----------------------------------------------------

# Imports.
using SpinAdaptedSecondQuantization


# Partitioning of electronic indices.
orbitalA    = new_space(:orbitalA, "A", "pqrstuv") 
orbitalB    = new_space(:orbitalB, "B", "pqrstuv")
virOrbitalA = new_space(:virOrbitalA, "voA", "abcdefg");
occOrbitalA = new_space(:occOrbitalA, "ooA", "ijklmno");
virOrbitalB = new_space(:virOrbitalB, "voB", "abcdefg");
occOrbitalB = new_space(:occOrbitalB, "ooB", "ijklmno");


# Enforce splitting s.t. orbitalsA, orbitalsB ⊆ GeneralOrbital, orbitalsA ∪ orbitalsB = GeneralOrbital, orbitalsA ∩ orbitalsB = {} etc.
add_space_sum(orbitalA, orbitalB, GeneralOrbital)
add_space_sum(virOrbitalA, virOrbitalB, VirtualOrbital)
add_space_sum(occOrbitalA, occOrbitalB, OccupiedOrbital)
add_space_sum(occOrbitalA, virOrbitalA, orbitalA)
add_space_sum(occOrbitalB, virOrbitalB, orbitalB)


# Enfore splitting s.t. orbitalA ∩ OccupiedOrbital = occOrbitalA etc. 
add_space_intersection(orbitalA, OccupiedOrbital, occOrbitalA)
add_space_intersection(orbitalA, VirtualOrbital, virOrbitalA)
add_space_intersection(orbitalB, OccupiedOrbital, occOrbitalB)
add_space_intersection(orbitalB, VirtualOrbital, virOrbitalB)


# Color coding such that indices from different spaces may be distinguished from each another, red := A, cyan := B. 
set_color(orbitalA, :red) 
set_color(orbitalB, :cyan)
set_color(virOrbitalA, :red)
set_color(occOrbitalA, :red)
set_color(virOrbitalB, :cyan)
set_color(occOrbitalB, :cyan)


# Definiting operators.
# One-electron parts.
gA  = real_tensor("F", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3])
gB  = real_tensor("F", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3])


hA  = summation(gA * E(1,2) * constrain(1 => orbitalA, 2 => orbitalA), 1:2)
hB  = summation(gB * E(1,2) * constrain(1 => orbitalB, 2 => orbitalB), 1:2)
hAB = summation(psym_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalA, 2 => orbitalB), 1:2)       
hBA = summation(psym_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalB, 2 => orbitalA), 1:2)      


# two-electron parts.
# particle conserving terms.
gA = 1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalA, 3 => orbitalA, 4 => orbitalA), 1:4)
gB = 1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalB, 3 => orbitalB, 4 => orbitalB), 1:4)
gBBAA = summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalB, 3 => orbitalA, 4 => orbitalA), 1:4)
gABBA = summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalB, 4 => orbitalA), 1:4)


# particle breaking terms.
gBAAA = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalA, 3 => orbitalA, 4 => orbitalA), 1:4))
gABAA = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalA, 4 => orbitalA), 1:4))
gABAB = simplify(1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalA, 4 => orbitalB), 1:4))
gBABA = simplify(1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalA, 3 => orbitalB, 4 => orbitalA), 1:4))
gABBB = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalB, 4 => orbitalB), 1:4))
gBABB = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalA, 3 => orbitalB, 4 => orbitalB), 1:4))


# Hamiltonian:
H_A  = hA + gA
H_B  = hB + gB
H_ABpc = gABBA + gBBAA


H_ABpb = hAB + hBA + gBAAA + gABAA + gABAB + gBABA + gABBB + gBABB
#H_ABpb = gBAAA + gABAA + gABAB + gBABA + gABBB + gBABB


Hpc  = H_A + H_B + H_ABpc
#Hpc  = H_ABpc

Hpb = H_ABpb
H   = Hpc + Hpb


# Excitation operators
T_A = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalA, 2 => occOrbitalA), 1:2)   
T_B = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalB, 2 => occOrbitalB), 1:2)
T_1 = T_A + T_B                                                                                                             # T2 = 0 to simplify derivations. 


# Functions for single-excitation matrix elements:
function nonest()
    Eia = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Ebj = E(3,4) * constrain(3 => virOrbitalB, 4 => occOrbitalA)                 


    # Verification of the first term (single commutator term)
    println("-------TERM 1--------")
    ket = commutator(Hpc,  Ebj)
    bra = 1//2 * Eia
    HFexpect = simplify_heavy(hf_expectation_value(bra * ket))

    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L")) 
    println("<HF Eia [Hpc, Ebj] HF> = $(HFexpect)")
    println()


    # Verification of sub terms (splittings A, B and AB) for the single commutator term
    bra = 1//2 * Eia
    ket = commutator(H_A, Ebj)
    println("1/2 * <HF Eia [H_A, Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")
    println()


    ket = commutator(H_B, Ebj)  
    println("1/2 * <HF Eia [H_B, Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")
    println()

    ket = commutator(H_ABpc, Ebj)
    println("1/2 * <HF Eia [H_ABpc, Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")
    println()
end


function onenest()
    Eia  = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Ebj  = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 


    # Verification of the second term (singly nested commutator term)
    println("-------TERM 2 -------")
    bra = 1//2 * Eia
    ket = commutator(Hpc, T_1, Ebj) 
    HFexpect = simplify_heavy(hf_expectation_value(bra*ket))
    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L"))
    println("1/2 * <HF Eia [[H_ABpc, T], Ebj] HF> = $(simplify_heavy(HFexpect)))")
    println()


    # Verification of sub terms (splittings A, B and AB) for the singly nested commutator term
    ket = commutator(H_A, T_A, Ebj) 
    HFexpect = simplify_heavy(hf_expectation_value(bra*ket)) 
    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L"))
    println("1/2 * <HF Eia [[H_A, T_A], Ebj] HF> = $(simplify_heavy(HFexpect))")
    println()


    ket = commutator(H_B, T_B, Ebj)
    HFexpect = simplify_heavy(hf_expectation_value(bra*ket))
    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L"))
    println("1/2 * <HF Eia [[H_B, T_B], Ebj] HF> = $(simplify_heavy(HFexpect))")
    println()


    ket = commutator(H_ABpc, T_1, Ebj)
    HFexpect = simplify_heavy(hf_expectation_value(bra*ket))
    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L"))
    println("1/2 * <HF Eia [[H_ABpc, T], Ebj] HF> = $(simplify_heavy(HFexpect))")
    println()
end


function twonest()
    Eia = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)
    Ebj = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)


    # Verification of the third term (doubly nested commutator term)
    println("-------TERM 3 -------")
    bra = 1//2 * Eia 
    
    ket = commutator(H_A, T_A, T_A, Ebj)
    HFexpect = simplify_heavy(hf_expectation_value(bra*ket))
    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L"))
    println("1/2 <HF Eia [[[H_A, T_A], T_A], Ebj] HF> = $(simplify_heavy(HFexpect))")
    println()


    ket = commutator(H_B, T_B, T_B, Ebj)
    HFexpect = simplify_heavy(hf_expectation_value(bra*ket))
    HFexpect = look_for_tensor_replacements(HFexpect, make_exchange_transformer("g", "L"))
    println("1/2 <HF Eia [[[H_B, T_B], T_B], EBj] HF > = $(simplify_heavy(HFexpect))")
    println()
end



# Testing:
function testing() 
    println("TESTING")
    Ebj = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)
    epqrs = e(3,4,5,6) * constrain(3 => orbitalB, 4 => orbitalB, 5 => virOrbitalA, 6 => occOrbitalB)
    

    println("[epqrs, Eem] = $(commutator(Ebj, epqrs))")

    #bra = Eik
    #ket = Ecm
    #println("<HF Eik Ecm HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")                                   
end


nonest()
onenest()
twonest()
#testing()

