# -----------------------------------------------------
# File: a11_verification.jl
# Author: Ivan Nygaard
# Created: 28/10/2025
# Description: Code used to verify the derivation of
#              the first element of the CC-Jacobian
#              that was derived by hand. 
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

# Debugging
#F_pq = real_tensor("F", 1,2) + summation((-2 * rsym_tensor("g", 1,2,3,3) + rsym_tensor("g", 1,3,3,2)) * constrain(3 => OccupiedOrbital), [3])
#hA = summation(real_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalA, 2 => orbitalA), 1:2)
#hB = summation(real_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalB, 2 => orbitalB), 1:2)
#FAeff  = real_tensor("FAeff", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3]) + summation((2 * rsym_tensor("g", 1,2,3,3) - rsym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3])
#FBeff  = real_tensor("FBeff", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3]) + summation((2 * rsym_tensor("g", 1,2,3,3) - rsym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3])

gA  = real_tensor("F", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3])
gB  = real_tensor("F", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3])


hA  = summation(gA * E(1,2) * constrain(1 => orbitalA, 2 => orbitalA), 1:2)
hB  = summation(gB * E(1,2) * constrain(1 => orbitalB, 2 => orbitalB), 1:2)
#hAB = summation(F_pq * E(1,2) * constrain(1 => orbitalA, 2 => orbitalB), 1:2)                       # F_AB
#hBA = summation(F_pq * E(1,2) * constrain(1 => orbitalB, 2 => orbitalA), 1:2)                       # F_AB
hAB = summation(psym_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalA, 2 => orbitalB), 1:2)       # F_AB
hBA = summation(psym_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalB, 2 => orbitalA), 1:2)       # F_BA


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
H_ABpc = gBBAA + gABBA


H_ABpb = hAB + hBA + gBAAA + gABAA + gABAB + gBABA + gABBB + gBABB
#H_ABpb = gBAAA + gABAA + gABAB + gBABA + gABBB + gBABB


Hpc  = H_A + H_B + H_ABpc
#Hpc  = H_ABpc:wq

Hpb = H_ABpb
H = Hpc + Hpb


# Excitation operators
T_A = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalA, 2 => occOrbitalA), 1:2)    # should be 1/4 for T2 but in overleaf uses 1/2? 
T_B = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalB, 2 => occOrbitalB), 1:2)
T_1 = T_A + T_B


# Functions for single-excitation matrix elements:
function nonest()
    Eia = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Ebj = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 

    # Verification of the first term (single commutator term)
    println("-------TERM 1--------")
    ket = commutator(Hpc,  Ebj)
    bra = 1//2 * Eia
    println("1/2 * <HF Eia [H(pc), Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")
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



    #braket = look_for_tensor_replacements(braket, make_exchange_transformer("g", "L")) 
    #braket = look_for_tensor_replacements(braket, make_exchange_transformer("t", "u"))

    #disable_external_index_translation()

    #braket = simplify_heavy(braket)
end


function onenest()
    Eia  = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Ebj  = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 


    # Verification of the second term (singly nested commutator term)
    println("-------TERM 2 -------")
    bra = 1//2 * Eia
    ket = commutator(Hpc, T_1, Ebj) 
    println("1/2 * <HF Eia [[H_ABpc, T], Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")


    # Verification of sub terms (splittings A, B and AB) for the singly nested commutator term
    ket = commutator(H_A, T_A, Ebj) 
    println("1/2 * <HF Eia [[H_A, T_A], Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")


    ket = commutator(H_B, T_B, Ebj)
    println("1/2 * <HF Eia [[H_B, T_B], Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")


    ket = commutator(H_ABpc, T_1, Ebj)
    println("1/2 * <HF Eia [[H_ABpc, T], Ebj] HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")
end




# Testing:
function testing()
    Eba = E(1,2) * constrain(1 => virOrbitalA, 2 => virOrbitalA)   
    Ers = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalA)                 
    Eai = E(5,6) * constrain(5 => occOrbitalB, 6 => virOrbitalA)

    println("TESTING")
    bra = Eba
    ket = Ers
    println("<HF E_baE_rs  HF> = $(simplify_heavy(hf_expectation_value(bra * ket)))")                                   
end


#nonest()
#onenest()
testing()

