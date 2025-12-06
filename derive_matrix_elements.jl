# -----------------------------------------------------
# File: derive_matrix_elements.jl
# Author: Ivan Nygaard
# Created: 05/11/2025
# Description: Code used to derive the CCSD-EOM
#              matrix (A, eta, omega). 
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
#gAeff  = real_tensor("FAeff", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3]) + summation((2 * psym_tensor("g", 1,2,3,3) - psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3])
#gBeff  = real_tensor("FBeff", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3]) + summation((2 * psym_tensor("g", 1,2,3,3) - psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3])


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
#Hpc  = H_ABpc
Hpb = H_ABpb
H = Hpc + Hpb


# Excitation operators
T_A = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalA, 2 => occOrbitalA), 1:2)    # should be 1/4 for T2 but in overleaf uses 1/2? 
T_B = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalB, 2 => occOrbitalB), 1:2)
T_1 = T_A + T_B


# Functions for single-excitation matrix elements:
function A11()
    Eia = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Ebj = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 

    bra = 1//2 * act_eT_on_bra(Eia, -T_1) * commutator(H, Ebj)
    braket = simplify_heavy(act_eT_on_bra(bra, T_1, max_ops = 0))


    braket = look_for_tensor_replacements(braket, make_exchange_transformer("g", "L")) 
    braket = look_for_tensor_replacements(braket, make_exchange_transformer("t", "u"))

    #disable_external_index_translation()

    braket = simplify_heavy(braket)
end


function A12()
    Eia = E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Ebj = E(3,4) * constrain(3 => virOrbitalB, 4 => occOrbitalA)                 

    bra = 1//2 * act_eT_on_bra(Eia, -T_1) * commutator(H, Ebj)
    braket = simplify_heavy(act_eT_on_bra(bra, T_1, max_ops = 0))
     
    
    braket = look_for_tensor_replacements(braket, make_exchange_transformer("g", "L")) 
    braket = simplify_heavy(braket)
end


function A21()
    Eia = E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    Ebj = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 

    bra = 1//2 * act_eT_on_bra(Eia, -T_1) * commutator(H, Ebj)
    braket = simplify_heavy(act_eT_on_bra(bra, T_1, max_ops = 0))

    braket = look_for_tensor_replacements(braket, make_exchange_transformer("g", "L")) 
end


function A22()
    Eia = E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    Ebj = E(3,4) * constrain(3 => virOrbitalB, 4 => occOrbitalA)                 

    bra = 1//2 * act_eT_on_bra(Eia, -T_1) * commutator(H, Ebj)
    braket = simplify_heavy(act_eT_on_bra(bra, T_1, max_ops = 0))
    braket = look_for_tensor_replacements(braket, make_exchange_transformer("g", "L")) 
    braket = simplify_heavy(braket)
end   



# Printing
println("A_11 = $(A11())")
println()
println("A_12 = $(A12())")
println()
println("A_21 = $(A21())")
println()
println("A_22 = $(A22())")
println()



# TO-DO: 
# 1. Define double excitation operators, note factor 1/4 in MEST! 
# 2. Check that single excitation terms remain the same when redefening T to be T1 + T2
# 3. Write functions for the double terms between doubly excited determinants in the overall matrix, also what about these eta terms? 
# 4. Can also check that all terms in the first column (except 11 = E0_CC) go to zero in the code too. 
