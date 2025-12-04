using SpinAdaptedSecondQuantization

h = summation(rsym_tensor("h", 1,2) * E(1,2) * electron(1,2), 1:2) 
g = simplify(1//2 * summation(rsym_tensor("g", 1:4...) * e(1:4...) * electron(1:4...), 1:4))   # Contracts delta (reduces number of summation indices)

H = h+g

E_hf = simplify_heavy(hf_expectation_value(H))						                            # Heavy also looks for integral symmetries to simplify further (in additions to killing dummies with kdeltas)

println(E_hf)
