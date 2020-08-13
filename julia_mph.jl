
#####################################################################################################################
# Maximum-peak height (MPH) Sentinel-3 OLCI Julia example implementation - Mortimer Werther (2020)                  #
# The code is based on the publication by M. Matthews and D. Odermatt 2015:                                         #
# "Improved algorithm for routine monitoring of cyanobacteria and eutrophication in inland and near-coastal waters" #
# https://doi.org/10.1016/j.rse.2014.10.010                                                                         #
#####################################################################################################################

# Splitting the calcuation of MPH CHL into two functions:
# 1. Required MPH run procedure parameters such as the peaks and positions
# 2. MPH run procedure using the parameters from 1.

function calc_mph_paras(rrs_620, rrs_665, rrs_681, rrs_709, rrs_753, rrs_885)

    """
    1. Calculation of required MPH parameters
    Requires satellite-derived BRR or in situ reflectance at 620nm, 665nm, 681nm, 709nm, 753nm and 885nm (OLCI)
    """

    # 1. calculate rmax_0, lambda_rmax_0, OLCI bands 8-9

    rmax_0 = max(rrs_681, rrs_709)

    if rrs_681 >= rrs_709
        lambda_rmax_0 = 681
    else
        lambda_rmax_0 = 709
    end

    # 2. calculate rmax_1, lambda_rmax_1, OLCI bands 8-10
    rmax_1 = max(rrs_681, rrs_709, rrs_753)

    # short-circuit evaluation using &&, can also use brackets around the arguments to use the single & operator
    if rrs_681 >= rrs_709 && rrs_681 >= rrs_753
        lambda_rmax_1 = 681
    elseif rrs_709 >= rrs_681 && rrs_709 >= rrs_753
        lambda_rmax_1 = 709
    else
        lambda_rmax_1 = 753
    end

    # 3. NDVI
    ndvi = (rrs_885 - rrs_665) / (rrs_885 + rrs_665)

    # 3. SIPF, SICF, BAIR
    sipf = rrs_665 - rrs_620 - ((rrs_681 - rrs_620) * (665 - 620) / (681 - 620))
    sicf = rrs_681 - rrs_665 - ((rrs_709 - rrs_665) * (681 - 665) / (709 - 665))
    bair = rrs_709 - rrs_665 - ((rrs_885 - rrs_665) * (709 - 665) / (885 - 665))

    # 4. MPH_0, MPH_1
    mph_0 = rmax_0 - rrs_665 - ((rrs_885 - rrs_665) * (lambda_rmax_0 - 665) / (885 - 665))
    mph_1 = rmax_1 - rrs_665 - ((rrs_885 - rrs_665) * (lambda_rmax_1 - 665) / (885 - 665))

    return mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0

end

# Definiton of MPH run procedure
function calc_mph_chl(mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0, mph_floatthres, mph_cyanomax)

    """
    Calculation of MPH CHL
    Requires previously calculated parameters AND mph_floatthres and mph_cyanomax (high int64 values, can be derived from paper)
    Feel free to uncomment or delete print statements.
    """

    if lambda_rmax_1 == 753
        println("Right side of if-condition.")
        if (mph_1 >= 0.02) | (ndvi >= 0.2)
            float_flag = 1
            adj_flag = 0
            # SICF < 0 and SIPF > 0
            if (sicf < 0) & (sipf > 0)
                cyano_flag=1
                println("Flag: floating cyanobacteria true")
                chl_mph = 22.44 * exp(35.79 * mph_1)
                # string concatenation is done with "*" in Julia - odd? consider hello^3 = hello * hello * hello
                if chl_mph > mph_floatthres
                    float_flag=1
                    println("Floating cyanobacteria")
                end
            # SICF >=0 or SIPF <=0
            elseif (sicf >= 0) | (sipf <= 0)
                cyano_flag = 0
                # equivalent to NaN in Python, NA in R etc.
                chl_mph = missing
                print("Floating vegetation")
            end
        # Continuation right side
        elseif (mph_1 < 0.02) & (ndvi < 0.2)
            print("continuation right side")
            float_flag = 0
            adj_flag = 1
            print("Flag: adjacent true")
            cyano_flag = 0
            print("Immersed eukaryotes")

            chl_mph = 5.24 * 10 ^ 9 * mph_0 ^ 4 - 1.95 * 10 ^ 8 * mph_0 ^ 3 + 2.46 * 10 ^ 6 * mph_0 ^ 2 + 4.02 * 10 ^ 3 * mph_0 + 1.97
        end
    else
        print("Left side of if-condition.")
        float_flag = 0
        adj_flag = 0

        # Left side of 2nd if-condition
        if (sicf >= 0) | (sipf <= 0) | (bair <= 0.002)
                print("Left 2nd if-condition")
                cyano_flag=0
                print("Immersed eukaryotes")
                chl_mph = 5.24 * 10 ^ 9 * mph_0 ^ 4 - 1.95 * 10 ^ 8 * mph_0 ^ 3 + 2.46 * 10 ^ 6 * mph_0 ^ 2 + 4.02 * 10 ^ 3 * mph_0 + 1.97

        # Right side of 2nd if-condition
        elseif (sicf <= 0) & (sipf > 0) & (bair > 0.002)
            print("Right 2nd if-condition")
            cyano_flag = 1
            print("Flag: cyanobacteria true")
            chl_mph = 22.44 * exp(35.79 * mph_1)
            if chl_mph > mph_floatthres
                    float_flag = 1
                    print("Floating cyanobacteria")
                    if chl_mph > mph_cyanomax
                        chl_mph = mph_cyanomax
                        print("MPH chl maximum reached.")
                    end
            end
        end
    end
    return chl_mph, adj_flag, cyano_flag, float_flag
end


# Sample Rrs (Sentinel - 3 OLCI bands)
rrs_620 = 0.0014
rrs_665 = 0.0020
rrs_681 = 0.0042
rrs_709 = 0.0025
rrs_753 = 0.0010
rrs_885 = 0.0003

# Call parameter function
mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0 = calc_mph_paras(rrs_620, rrs_665, rrs_681, rrs_709, rrs_753, rrs_885)

# Call MPH run procedure
chl_mph, adj_flag, cyano_flag, float_flag = calc_mph_chl(mph_0, mph_1, sipf, sicf, bair, ndvi, rmax_1, lambda_rmax_1, lambda_rmax_0, rmax_0, 200, 350)

println("MPH CHL is: " * string(chl_mph))
println("Flags - Adj flag: " * string(adj_flag)* ", "* "Cyano Flag: " * string(cyano_flag) * ", "* "Float flag: " * string(floats_flag))
