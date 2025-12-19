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


# Debugging
#F_pq = real_tensor("F", 1,2) + summation((-2 * rsym_tensor("g", 1,2,3,3) + rsym_tensor("g", 1,3,3,2)) * constrain(3 => OccupiedOrbital), [3])
#hA = summation(real_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalA, 2 => orbitalA), 1:2)
#hB = summation(real_tensor("h", 1,2) * E(1,2) * constrain(1 => orbitalB, 2 => orbitalB), 1:2)
#gAeff  = real_tensor("FAeff", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3]) + summation((2 * psym_tensor("g", 1,2,3,3) - psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3])
#gBeff  = real_tensor("FBeff", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3]) + summation((2 * psym_tensor("g", 1,2,3,3) - psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3])
#disable_external_index_translation()
#gA  = real_tensor("F", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalA), [3])
#gB  = real_tensor("F", 1,2) + summation((-2 * psym_tensor("g", 1,2,3,3) + psym_tensor("g", 1,3,3,2)) * constrain(3 => occOrbitalB), [3])
#hAB = summation(F_pq * E(1,2) * constrain(1 => orbitalA, 2 => orbitalB), 1:2)                       
#hBA = summation(F_pq * E(1,2) * constrain(1 => orbitalB, 2 => orbitalA), 1:2)                       


# Defining operatos
# One-electron parts.
gAeff  = real_tensor("Feff", 1,2)  - summation(psym_tensor("L", 1,2,3,3) * constrain(3 => occOrbitalA), [3]) - summation(psym_tensor("L", 1,2,4,4) * constrain(4 => occOrbitalB), [4])
gBeff  = real_tensor("Feff", 1,2)  - summation(psym_tensor("L", 1,2,3,3) * constrain(3 => occOrbitalB), [3]) - summation(psym_tensor("L", 1,2,4,4) * constrain(4 => occOrbitalA), [4])
gABeff = real_tensor("Feff", 1,2)  - summation(psym_tensor("L", 1,2,3,3) * constrain(3 => OccupiedOrbital), [3])                                                                                # Intuisjon om occ here, men why? 
gBAeff = gABeff

hA  = summation(gAeff * E(1,2) * constrain(1 => orbitalA, 2 => orbitalA), 1:2)
hB  = summation(gBeff * E(1,2) * constrain(1 => orbitalB, 2 => orbitalB), 1:2)
hAB = summation(gABeff * E(1,2) * constrain(1 => orbitalA, 2 => orbitalB), 1:2)       
hBA = summation(gBAeff * E(1,2) * constrain(1 => orbitalB, 2 => orbitalA), 1:2)       


# Two-electron parts.
# Particle conserving terms. 
gA = 1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalA, 3 => orbitalA, 4 => orbitalA), 1:4)
gB = 1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalB, 3 => orbitalB, 4 => orbitalB), 1:4)
gBBAA = summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalB, 3 => orbitalA, 4 => orbitalA), 1:4) 
gABBA = summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalB, 4 => orbitalA), 1:4)


# Particle breaking terms.
gBAAA = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalA, 3 => orbitalA, 4 => orbitalA), 1:4))
gABAA = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalA, 4 => orbitalA), 1:4))
gABAB = simplify(1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalA, 4 => orbitalB), 1:4))
gBABA = simplify(1//2 * summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalA, 3 => orbitalB, 4 => orbitalA), 1:4))
gABBB = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalA, 2 => orbitalB, 3 => orbitalB, 4 => orbitalB), 1:4))
gBABB = simplify(summation(psym_tensor("g", 1:4...) * e(1:4...) * constrain(1 => orbitalB, 2 => orbitalA, 3 => orbitalB, 4 => orbitalB), 1:4))


# Hamiltonian.
H_A  = hA + gA
H_B  = hB + gB
H_ABpc = gBBAA + gABBA

H_ABpb = hAB + hBA + gBAAA + gABAA + gABAB + gBABA + gABBB + gBABB

Hpc = H_A + H_B + H_ABpc
Hpb = H_ABpb
H   = Hpc + Hpb


# Excitation operators
T1_A  = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalA, 2 => occOrbitalA), 1:2)    
T1_B  = summation(real_tensor("t", 1,2) * E(1,2) * constrain(1 => virOrbitalB, 2 => occOrbitalB), 1:2)
T1_AB = 0                                                                                                                                                                # modified version, satisfies [T,n_A] = 0 = [T,n_B]
T2_A  = 1//2 * summation(real_tensor("t", 1:4...) * E(1,2) * E(3,4) * constrain(1 => virOrbitalA, 2 => occOrbitalA, 3 => virOrbitalA, 4 => occOrbitalA), 1:4)
T2_B  = 1//2 * summation(real_tensor("t", 1:4...) * E(1,2) * E(3,4) * constrain(1 => virOrbitalB, 2 => occOrbitalB, 3 => virOrbitalB, 4 => occOrbitalB), 1:4)
T2_AB = 1//2 * (summation(real_tensor("t", 1:4...) * E(1,2) * E(3,4) * constrain(1 => virOrbitalA, 2 => occOrbitalB, 3 => virOrbitalB, 4 => occOrbitalA), 1:4) +         # modified version, satisfies [T,n_A] = 0 = [T,n_B]
                summation(real_tensor("t", 1:4...) * E(1,2) * E(3,4) * constrain(1 => virOrbitalB, 2 => occOrbitalB, 3 => virOrbitalA, 4 => occOrbitalA), 1:4))


T_1 = T1_A + T1_B + T1_AB
T_2 = T2_A + T2_B + T2_AB
T   = T_2 # + T_1                                                                                                                                                        # T1 is commented out because we are operating with T1 transformed integrals (reduced symmetry) 


# Functions used to construct the matrix elements of the shifted CCSD-EOM Hamiltonian based on equations (), () and (). 
# i, j, k, l := occ, a,b,c,d = vir. In the name of the function that produces the defining equation for the relevant matrix element, bar and non-bar indices
# are explicitly stated in the name as the index followed by bar e.g. ibar for occupied orbital localized to B. 
function generate_pb_resolution_of_the_identity_states()
    E1 = E(9,10)  *            constrain(9  =>  virOrbitalA, 10 => occOrbitalB)
    E2 = E(11,12) *            constrain(11 =>  virOrbitalB, 12 => occOrbitalA) 
    E3 = E(13,14) * E(15,16) * constrain(13 =>  virOrbitalB, 14 => occOrbitalA, 15 => virOrbitalA, 16 => occOrbitalA) 
    E4 = E(17,18) * E(19,20) * constrain(17 =>  virOrbitalA, 18 => occOrbitalB, 19 => virOrbitalA, 20 => occOrbitalA)
    E5 = E(21,22) * E(23,24) * constrain(21 =>  virOrbitalA, 22 => occOrbitalB, 23 => virOrbitalA, 24 => occOrbitalB)
    E6 = E(25,26) * E(27,28) * constrain(25 =>  virOrbitalB, 26 => occOrbitalA, 27 => virOrbitalB, 28 => occOrbitalA)
    E7 = E(29,30) * E(31,32) * constrain(29 =>  virOrbitalB, 30 => occOrbitalB, 31 => virOrbitalA, 32 => occOrbitalB)
    E8 = E(33,34) * E(35,36) * constrain(33 =>  virOrbitalB, 34 => occOrbitalB, 35 => virOrbitalB, 36 => occOrbitalA)


    E3n = E(13,16) * E(15,14) * constrain(13 => virOrbitalB, 14 => occOrbitalA, 15 => virOrbitalA, 16 => occOrbitalA)  # Needed for normalization <aibjI = 1/3 <aibjI + 1/6 <ajbiI
    E4n = E(17,20) * E(19,18) * constrain(17 => virOrbitalA, 18 => occOrbitalB, 19 => virOrbitalA, 20 => occOrbitalA)
    E5n = E(21,24) * E(23,22) * constrain(21 => virOrbitalA, 22 => occOrbitalB, 23 => virOrbitalA, 24 => occOrbitalB)
    E6n = E(25,28) * E(27,26) * constrain(25 => virOrbitalB, 26 => occOrbitalA, 27 => virOrbitalB, 28 => occOrbitalA)
    E7n = E(29,32) * E(31,30) * constrain(29 => virOrbitalB, 30 => occOrbitalB, 31 => virOrbitalA, 32 => occOrbitalB)
    E8n = E(33,36) * E(35,34) * constrain(33 => virOrbitalB, 34 => occOrbitalB, 35 => virOrbitalB, 36 => occOrbitalA)


    kets = [E1, E2, E3, E4, E5, E6, E7, E8]
    bras = [1//2 * E1', 1//2 * E2', 1//3 * E3' + 1//6 * E3n', 1//3 * E4' + 1//6 * E4n', 1//3 * E5' + 1//6 * E5n', 1//3 * E6' + 1//6 * E6n', 1//3 * E7' + 1//6 * E7n', 1//3 * E8' + 1//6 * E8n']
    

    [bras, kets]
end


# Functions for the single-single matrix elements of the level-shifted CCSD-Jacobian:
function A_a_ibar_c_kbar()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Eck =        E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 

    bra = act_eT_on_bra(Eia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))

    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))

    res
end


function A_abar_i_c_kbar()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    Eck = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 

    bra = act_eT_on_bra(Eia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))

    res
end


function A_a_ibar_cbar_k()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    Eck = E(3,4) * constrain(3 => virOrbitalB, 4 => occOrbitalA)                 

    bra = act_eT_on_bra(Eia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))


    res
end


function A_abar_i_cbar_k()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    Eck = E(3,4) * constrain(3 => virOrbitalB, 4 => occOrbitalA)                 

    bra = act_eT_on_bra(Eia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))


    res
end


# Functions for the single-double matrix elements of the level-shifted CCSD-Jacobian. 
function A_a_ibar_c_kbar_d_l()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalA, 4 => occOrbitalB, 5 => virOrbitalA, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_cbar_kbar_d_lbar() 
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalB, 5 => virOrbitalA, 6 => occOrbitalB)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_cbar_k_d_l()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalA, 5 => virOrbitalA, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_cbar_kbar_dbar_l() 
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalB, 5 => virOrbitalB, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_c_kbar_d_lbar() 
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalA, 4 => occOrbitalB, 5 => virOrbitalA, 6 => occOrbitalB)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_cbar_k_dbar_l()  
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalA, 5 => virOrbitalB, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_c_kbar_d_l()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalA, 4 => occOrbitalB, 5 => virOrbitalA, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_cbar_kbar_d_lbar() 
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalB, 5 => virOrbitalA, 6 => occOrbitalB)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_cbar_k_d_l()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalA, 5 => virOrbitalA, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_cbar_kbar_dbar_l()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalB, 5 => virOrbitalB, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_c_kbar_d_lbar()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalA, 4 => occOrbitalB, 5 => virOrbitalA, 6 => occOrbitalB)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_cbar_k_dbar_l()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, 2 => virOrbitalB)   
    EckEdl = E(3,4) * E(5,6) * constrain(3 => virOrbitalB, 4 => occOrbitalA, 5 => virOrbitalB, 6 => occOrbitalA)                 


    bra = act_eT_on_bra(Eia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(Eia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(Eia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(Eia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


# Functions for the double-single matrix elements of the level-shifted CCSD-Jacobian. 
function A_a_ibar_b_j_c_kbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    Eck = E(5,6) * constrain(5 => virOrbitalA, 6 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_c_kbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalA, 6 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_c_kbar()  
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalA, 6 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_c_kbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalA, 6 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_c_kbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    Eck = E(5,6) * constrain(5 => virOrbitalA, 6 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_c_kbar()  
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalA, 6 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_j_cbar_k() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    Eck = E(5,6) * constrain(5 => virOrbitalB, 6 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_cbar_k() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalB, 6 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_cbar_k()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalB, 6 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_cbar_k()  
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalB, 6 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_cbar_k() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    Eck = E(5,6) * constrain(5 => virOrbitalB, 6 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_cbar_k()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    Eck = E(5,6) * constrain(5 => virOrbitalB, 6 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, Eck)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * Eck * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * Eck * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * Eck * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


# Functions for the double-double matrix elements of the level-shifted CCSD-Jacobian. 
function A_a_ibar_b_j_c_kbar_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_j_cbar_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_j_cbar_k_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_j_cbar_kbar_dbar_l()  
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_j_c_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_j_cbar_k_dbar_l() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_c_kbar_d_l() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_cbar_kbar_d_lbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_cbar_k_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_cbar_kbar_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_c_kbar_d_lbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_b_jbar_cbar_k_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_c_kbar_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_cbar_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_cbar_k_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_cbar_kbar_dbar_l() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_c_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_b_j_cbar_k_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_c_kbar_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_cbar_kbar_d_lbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_cbar_k_d_l() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_cbar_kbar_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_c_kbar_d_lbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_ibar_bbar_j_cbar_k_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_c_kbar_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_cbar_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_cbar_k_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_cbar_kbar_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_c_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_a_ibar_b_jbar_cbar_k_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_c_kbar_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_cbar_kbar_d_lbar() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_cbar_k_d_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_cbar_kbar_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_c_kbar_d_lbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function A_abar_i_bbar_j_cbar_k_dbar_l()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(EjbEia, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    bras = generate_pb_resolution_of_the_identity_states()[1]
    kets = generate_pb_resolution_of_the_identity_states()[2]

    braket2 = simplify_heavy(hf_expectation_value(EjbEia * EckEdl * kets[1]) * act_eT_on_bra((act_eT_on_bra(bras[1], -T) * H), T, max_ops = 0) +    
                             hf_expectation_value(EjbEia * EckEdl * kets[2]) * act_eT_on_bra((act_eT_on_bra(bras[2], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[3]) * act_eT_on_bra((act_eT_on_bra(bras[3], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[4]) * act_eT_on_bra((act_eT_on_bra(bras[4], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[5]) * act_eT_on_bra((act_eT_on_bra(bras[5], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[6]) * act_eT_on_bra((act_eT_on_bra(bras[6], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[7]) * act_eT_on_bra((act_eT_on_bra(bras[7], -T) * H), T, max_ops = 0) +
                             hf_expectation_value(EjbEia * EckEdl * kets[8]) * act_eT_on_bra((act_eT_on_bra(bras[8], -T) * H), T, max_ops = 0))


    res = braket1 + braket2
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


# Omega terms
function O_a_ibar()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalB, 2 => virOrbitalA)   

    bra = act_eT_on_bra(Eia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_a_ibar_b_j() 
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   
    
    bra = act_eT_on_bra(EjbEia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_abar_ibar_b_jbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalB)   
    bra = act_eT_on_bra(EjbEia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_abar_i()
    Eia = 1//2 * E(1,2) * constrain(1 => occOrbitalA, virOrbitalB)

    bra = act_eT_on_bra(Eia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_abar_i_b_j()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalA, 3 => occOrbitalA, 4 => virOrbitalB)   

    bra = act_eT_on_bra(EjbEia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_abar_ibar_bbar_j()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalB, 4 => virOrbitalB)   

    bra = act_eT_on_bra(EjbEia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_a_ibar_b_jbar()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalB, 2 => virOrbitalA, 3 => occOrbitalB, 4 => virOrbitalA)   

    bra = act_eT_on_bra(EjbEia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function O_abar_i_bbar_j()
    EjbEia = (1//3 * E(1,2) * E(3,4) + 1//6 * E(3,2) * E(1,4)) * constrain(1 => occOrbitalA, 2 => virOrbitalB, 3 => occOrbitalA, 4 => virOrbitalB)   
    bra = act_eT_on_bra(EjbEia, -T) * H
    braket1 = simplify_heavy(act_eT_on_bra(bra, T, max_ops = 0))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


# Eta terms, NB: Unfortunately, these do not work yet as I have not been able to tell SASQ that i want the bra to just be the HF-state. I tried id = 1,
# but act_eT_on_bra requires and expr, so i tread E_ii unsure if this is the right way to do this. 
function N_c_kbar()
    id = E(1,1) * constrain(1 => OccupiedOrbital, 2 => OccupiedOrbital)
    Eck = E(3,4) * constrain(3 => virOrbitalA, 4 => occOrbitalB)                 

    bra = act_eT_on_bra(id, -T) * commutator(H, Eck)                                                         # <HF E_ai = 0 so the <HF taumu^dagger tau_nu H^ST HF> = 0
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_c_kbar_d_l()
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(1, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_cbar_kbar_d_lbar() 
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)

    bra = act_eT_on_bra(1, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_cbar_k()
    Eck = E(3,4) * constrain(3 => virOrbitalB, 4 => occOrbitalA)

    bra = act_eT_on_bra(1, -T) * commutator(H, Eck)                                                        
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1 
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_cbar_k_d_l() 
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalA, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(1, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_cbar_kbar_dbar_l()
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalB, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(1, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_c_kbar_d_lbar()
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalA, 6 => occOrbitalB, 7 => virOrbitalA, 8 => occOrbitalB)                 

    bra = act_eT_on_bra(1, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


function N_cbar_k_dbar_l() 
    EckEdl = E(5,6) * E(7,8) * constrain(5 => virOrbitalB, 6 => occOrbitalA, 7 => virOrbitalB, 8 => occOrbitalA)                 

    bra = act_eT_on_bra(1, -T) * commutator(H, EckEdl)
    braket1 = simplify_heavy(act_eT_on_bra(bra, T))

    res = braket1
    res = look_for_tensor_replacements(res, make_exchange_transformer("g", "L"))
    res = look_for_tensor_replacements(res, make_exchange_transformer("t", "u"))

    res
end


# Display terms
# Single-single elements CCSD-Jacobian, 4
println("-------SINGLE-SINGLE ELEMENTS CCSD-JACOBIAN-------")
println("Aaĩ,ck̃ = $(A_a_ibar_c_kbar())")
println()
println("Aãi,ck̃ = $(A_abar_i_c_kbar())")
println()
println("Aaĩ,c̃k = $(A_a_ibar_cbar_k())")
println()
println("Aãi,c̃k = $(A_abar_i_cbar_k())")
println()


# Single-double elements CCSD-Jacobian, 12
println("-------SINGLE-DOUBLE ELEMENTS CCSD-JACOBIAN-------")
println("Aaĩ,ck̃dl = $(A_a_ibar_c_kbar_d_l())")
println()
println("Aaĩ,c̃k̃dl̃ = $(A_a_ibar_cbar_kbar_d_lbar())")
println()
#println("Aaĩ,c̃kdl = $(A_a_ibar_cbar_k_d_l())")
#println("Aaĩ,c̃k̃d̃l = $(A_a_ibar_cbar_kbar_dbar_l())")
#println("Aaĩ,ck̃dl̃ = $(A_a_ibar_c_kbar_d_lbar())")
#println("Aaĩ,c̃kd̃l = $(A_a_ibar_cbar_k_dbar_l)")
#println("Aãi,ck̃dl = $(A_abar_i_c_kbar_d_l())")
#println("Aãi,c̃k̃dl̃ = $(A_abar_i_cbar_kbar_d_lbar())")
#println("Aãi,c̃kdl = $(A_abar_i_cbar_k_d_l())")
#println("Aãi,c̃k̃d̃l = $(A_abar_i_cbar_kbar_dbar_l())")
#println("Aãi,ck̃dl̃ = $(A_abar_i_c_kbar_d_lbar())")
#println("Aãi,c̃kd̃l = $(A_abar_i_cbar_k_dbar_l())")


# Double-single elements CCSD-Jacobian, 12
#println("-------DOUBLE-SINGLE ELEMENTS CCSD-JACOBIAN-------")
#println("Aaĩbj,ck̃ = $(A_a_ibar_b_j_c_kbar())")
#println("Aãĩbj̃,ck̃ = $(A_abar_ibar_b_jbar_c_kbar())")
#println("Aãibj,ck̃ = $(A_abar_i_b_j_c_kbar())")
#println("Aãĩb̃j,ck̃ = $(A_abar_ibar_bbar_j_c_kbar())")
#println("Aaĩbj̃,ck̃ = $(A_a_ibar_b_jbar_c_kbar())")
#println("Aãib̃j,ck̃ = $(A_abar_i_bbar_j_c_kbar())")
#println("Aaĩbj,c̃k = $(A_a_ibar_b_j_cbar_k())")
#println("Aãĩbj̃,c̃k = $(A_abar_ibar_b_jbar_cbar_k())")
#println("Aãibj,c̃k = $(A_abar_i_b_j_cbar_k())")
#println("Aãĩb̃j,c̃k = $(A_abar_ibar_bbar_j_cbar_k())")
#println("Aaĩbj̃,c̃k = $(A_a_ibar_b_jbar_cbar_k())")
#println("Aãib̃j,c̃k = $(A_abar_i_bbar_j_cbar_k())")


# Double-double elements CCSD-Jacobian, 36
#println("-------DOUBLE-DOUBLE ELEMENTS CCSD-JACOBIAN-------")
#println("Aaĩbj,ck̃dl = $(A_a_ibar_b_j_c_kbar_d_l())")
#println()
#println("Aaĩbj,c̃k̃dl̃ = $(A_a_ibar_b_j_cbar_kbar_d_lbar())")
#println("Aaĩbj,c̃kdl = $(A_a_ibar_b_j_cbar_k_d_l())")
#println("Aaĩbj,c̃k̃d̃l = $(A_a_ibar_b_j_cbar_kbar_dbar_l())")
#println("Aaĩbj,ck̃dl̃ = $(A_a_ibar_b_j_c_kbar_d_lbar())")
#println("Aaĩbj,c̃kd̃l = $(A_a_ibar_b_j_cbar_k_dbar_l())")
#println("Aãĩbj̃,ck̃dl = $(A_abar_ibar_b_jbar_c_kbar_d_l())")
#println("Aãĩbj̃,c̃k̃dl̃ = $(A_abar_ibar_b_jbar_cbar_kbar_d_lbar())")
#println("Aãĩbj̃,c̃kdl = $(A_abar_ibar_b_jbar_cbar_k_d_l())")
#println("Aãĩbj̃,c̃k̃d̃l = $(A_abar_ibar_b_jbar_cbar_kbar_dbar_l())")
#println("Aãĩbj̃,ck̃dl̃ = $(A_abar_ibar_b_jbar_c_kbar_d_lbar())")
#println("Aãĩbj̃,c̃kd̃l = $(A_abar_ibar_b_jbar_cbar_k_dbar_l())")
#println("Aãibj,ck̃dl = $(A_abar_i_b_j_c_kbar_d_l())")
#println("Aãibj,c̃k̃dl̃ = $(A_abar_i_b_j_cbar_kbar_d_lbar())")
#println("Aãibj,c̃kdl = $(A_abar_i_b_j_cbar_k_d_l())")
#println("Aãibj,c̃k̃d̃l = $(A_abar_i_b_j_cbar_kbar_dbar_l())")
#println("Aãibj,ck̃dl̃ = $(A_abar_i_b_j_c_kbar_d_lbar())")
#println("Aãibj,c̃kd̃l = $(A_abar_i_b_j_cbar_k_dbar_l())")
#println("Aãĩb̃j,ck̃dl = $(A_abar_ibar_bbar_j_c_kbar_d_l())") 
#println("Aãĩb̃j,c̃k̃dl = $(A_abar_ibar_bbar_j_cbar_kbar_d_lbar())")
#println("Aãĩb̃j,c̃kdl = $(A_abar_ibar_bbar_j_cbar_k_d_l())")
#println("Aãĩb̃j,c̃k̃d̃l = $(A_abar_ibar_bbar_j_cbar_kbar_dbar_l())")
#println("Aãĩb̃j,ck̃dl̃ = $(A_abar_ibar_bbar_j_c_kbar_d_lbar())")
#println("Aãĩb̃j,c̃kd̃l = $(A_abar_ibar_bbar_j_cbar_k_dbar_l())")
#println("Aaĩbj̃,ck̃dl = $(A_a_ibar_b_jbar_c_kbar_d_l())")
#println("Aaĩbj̃,c̃k̃dl̃ = $(A_a_ibar_b_jbar_cbar_kbar_d_lbar())")
#println("Aaĩbj̃,c̃kdl = $(A_a_ibar_b_jbar_cbar_k_d_l())")
#println("Aaĩbj̃,c̃k̃d̃l = $(A_a_ibar_b_jbar_cbar_kbar_dbar_l())")
#println("Aaĩbj̃_ck̃dl̃ = $(A_a_ibar_b_jbar_c_kbar_d_lbar())")
#println("Aaĩbj̃,c̃kd̃l = $(A_a_ibar_b_jbar_cbar_k_dbar_l())")
#println("Aãib̃j,ck̃dl = $(A_abar_i_bbar_j_c_kbar_d_l())")
#println("Aãib̃j,c̃k̃dl̃ = $(A_abar_i_bbar_j_cbar_kbar_d_lbar())")
#println("Aãib̃j,c̃kdl = $(A_abar_i_bbar_j_cbar_k_d_l())")
#println("Aãib̃j,c̃k̃d̃l = $(A_abar_i_bbar_j_cbar_kbar_dbar_l())")
#println("Aãib̃j,ck̃dl̃ = $(A_abar_i_bbar_j_c_kbar_d_lbar())")
#println("Aãib̃j,c̃kd̃l = $(A_abar_i_bbar_j_cbar_k_dbar_l())")
#println()


# Omega, 7
#println("-------ELEMENTS OF OMEGA-VECTOR-------")
println("Oaĩ   = $(O_a_ibar())")
println()
println("Oaĩbj = $(O_a_ibar_b_j())")
println()
#println("Oãĩbj̃ = $(O_abar_ibar_b_jbar())")
#println("Oãi   = $(O_abar_i())")
#println("Oãibj = $(O_abar_i_b_j())")
#println("Oãĩb̃j = $(O_abar_ibar_bbar_j())")
#println("Oaĩbj̃ = $(O_a_ibar_b_jbar())")
#println("Oãib̃j = $(O_abar_i_bbar_j())")



# NB: These do unfortunately not work yet. 
# Eta, 7
#println("-------ELEMENTS OF ETA-VECTOR-------")
#println("Nck̃   = $(N_c_kbar())")
#println()
#println("Nck̃dl = $(N_c_kbar_d_l())")
#println("Nc̃k̃dl̃ = $(N_cbar_kbar_d_lbar())")
#println("Nc̃k   = $(N_cbar_k())")
#println("Nc̃kdl = $(N_cbar_k_d_l())")
#println("Nc̃k̃d̃l = $(N_cbar_kbar_dbar_l())")
#println("Nck̃dl̃ = $(N_c_kbar_d_lbar())")
#println("Nc̃kd̃l = $(N_cbar_k_dbar_l())")


# E0CC, 1
#println("-------E0_CC-------")
#println("E0CC $(E0CC())")
