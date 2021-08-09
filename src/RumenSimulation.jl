module RumenSimulation

using NLsolve
using Trapz

# Declare the initial concentration of ingested feed
C0 = 60 # g/L

# Declare a struct to hold the final results
struct RumenResults
    C
    k
    τ
end # struct

function cornell_model(BW, digestibility)
    # Calculate the volume of each digestive compartment (L)
    V = cornell_sizes(BW)

    # Calculate the intake on mass basis (kg/day)
    Q_m = (1/7) * V[1]

    # Calculate the intake on volume basis (L/day)
    Q = Q_m * 1e3 / C0

    # Calculate the mean residence time for each of the compartments
    τ = V ./ Q

    # Guess C and k based on "normal" literature values
    # See Dias et al. 2011
    C = [50.5 41.8 40.2 37.1 24.0]
    k = [0.054 0.054 0.054 0.054;
       0.056 0.056 0.056 0.056;
       0.083 0.083 0.083 0.083;
       0.166 0.166 0.166 0.166;
       0.123 0.123 0.123 0.123]
    k = k .* 24

    # Set up a one-parameter function to optimize
    optimfun(k) = cornell_outside_objective(C, k, τ, digestibility)

    # Calculate the best fit digestive rate parameters
    k_calc = nlsolve(optimfun, k, ftol=1e-3, show_trace=true)

    # Resolve for the best concentrations
    optimfun2(c) = cornell_inside_objective(c, k_calc.zero, τ)
    C_calc = nlsolve(optimfun2, C, ftol=1e-3, show_trace=true)

    # Return the results
    return RumenResults(C_calc.zero, k_calc.zero, τ)


end # function

function cornell_inside_objective(C, k, τ)
    # Set up a vector with the calculated residence times
    τ_calc = Array{Float64,2}(undef,1,5)

    # Calculate the expected residence time for each compartment
    for n in 1:5
        τ_calc[n] = cornell_design(n, k, C)
    end # for

    # Calculate the sum squared error
    errs = τ .- τ_calc
    # errs = sum(errs .^ 2)

    return errs
end # function

function cornell_outside_objective(C, k, τ, digestability)
    # Set up a one-parameter function to optimize
    optimfun(c) = cornell_inside_objective(c, k, τ)

    # Find the exit concentrations of each compartment
    C_calc = nlsolve(optimfun, C, ftol=1e-3)

    # Compare the fecal concentration to the given whole-
    # tract digestibility
    digestibility_calc = (C0 - C_calc.zero[end]) / C0
    err = digestibility_calc - digestability
    return err
end

function cornell_design(n, k, C)
    if n < 4
        return rumen_design(n, k, C)
    else
        return intestine_design(n, k, C)
    end # if
end # function

function cornell_sizes(BW)
    # Define the conversion from body weight to organ volume for each
    # compartment (L/kg)
    proportions = [6.95 4.266 4.582 0.5 4.6] ./ 100

    # Return the volumes based on those proportions
    return proportions .* BW
end

function rumen_design(n, k, C)
    # Extract the k parameters for better readability
    kd = k[n, 1]
    kr = k[n, 2]
    kp = k[n, 3]

    # Append the initial feed concentration to the concentration matrix
    C = [C0 C]

    # Remove passage if this is the first rumen compartment
    if n == 1
        kp = 0
    end

    # If this is under Michalis-Menten kinetics, the k vector will be larger
    if size(k, 2) == 4
        # We are in Michalis-Menten

        # Extract the dmax parameter
        dmax = k[n,4]

        # Here is the Michalis-Menten design equation
        τ_calc = (C[n] - C[n+1]) / (C[n+1] * (kr + kp) + ((dmax*C[n+1])/(kd + C[n+1])))

    else
        # We are in first-order regime

        # Here is the first-order design equation
        τ_calc = (C[n-2] - C[n-1]) / (C[n-1] * (kr + kp + kd))

    end #if

    return τ_calc
end # function

function intestine_design(n, k, C)
    # Extract the parameters
    kd = k[n,1]
    kp = k[n,3]
    dmax = k[n,4]

    # Find the upper limit of integration
    X2 = (C[n-1] - C[n]) / C[n-1]

    # Create a mesh of points from 0 to X2
    C_mesh = range(C[n-1], stop=X2, length=100)

    # Evaluate the integral expression at all points along the mesh
    τ_mesh = intestine_expression.(C_mesh, kd, dmax, kp)

    # Evaluate the integral
    τ_calc = trapz(C_mesh, τ_mesh)

    return τ_calc

end # function

function intestine_expression(C, kd, dmax, kp)
    return -1 / ((dmax*C)/(kd + C) + kp*C)
end # function

end # module
