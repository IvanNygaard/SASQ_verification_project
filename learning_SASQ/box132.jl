using SpinAdaptedSecondQuantization

h = summation((real_tensor("F", 1, 2) + 
	       summation((-2 * psym_tensor("g", 1, 2, 3, 3) + 
			       psym_tensor("g", 1, 3, 3, 2)) * occupied(3), [3])
			) * E(1, 2) * electron(1, 2), 1:2);

g = 1//2 * simplify(
           summation(psym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4)
       )

H = h + g									                                                # Every expression except the HF energy looks cleaner when expressed using Fock matrix elements(?).

Eai(a,i) = E(a,i) * virtual(a) * occupied(i)

H_HF       = simplify_heavy(act_on_ket(H))					                                # act_on_ket(A) := A┃HF>
H_E_HF     = simplify_heavy(act_on_ket(commutator(H, Eai(1,2))))                            # commutator(A,B) := [A,B]
H_E_E_HF   = simplify_heavy(act_on_ket(commutator(H, Eai(1,2), Eai(3,4))))	                # commutator(A,B,C) := [[A,B], C]
H_E_E_E_HF = simplify_heavy(act_on_ket(commutator(H, Eai(1,2), Eai(3,4), Eai(5,6))))        # commutator(A,B,C,D) := [[[A,B], C], D]

println("H = $H")
println()
println("H_HF = $H_HF")
println()

# Can simplify further by looking for L_pqrs = 2g_pqrs - g_psrq
H_E_HF = look_for_tensor_replacements(H_E_HF, make_exchange_transformer("g", "L"))

println("H_E_HF = $H_E_HF")
println()


# Can simplify further by using permutation operators (permutational symmetries stem from the Jacobi identity)
H_E_E_HF = look_for_tensor_replacements(H_E_E_HF, make_exchange_transformer("g", "L"))
r, ss, ns = desymmetrize(H_E_E_HF, make_permutation_mappings([(1, 2), (3, 4)]));

print("H_E_E_HF = $ss + $ns + $r")						      # symmetric terms + non-symmetric terms + P_ij^ab(...) terms 
