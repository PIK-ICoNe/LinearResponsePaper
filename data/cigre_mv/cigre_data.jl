using PowerDynamics: AbstractLine

function CIGRE_static(;close_switches)

    ###### knoten ######

    LoadP = [
        19839E3,
        0.0,
        501.7E3,
        431.65E3,
        727.5E3,
        548.05E3,
        76.5E3,
        586.85E3,
        573.75E3,
        543.3E3,
        329.8E3,
    ]

    GenP = [
        0.0,
        0.0,
        20E3,
        20E3,
        663E3,
        30E3,
        1500E3,
        30E3,
        552E3,
        254E3,
        10E3,
    ]

    cosϕ = [
        0.9737554,
        0.0,
        0.9231825,
        0.9700009,
        0.9699996,
        0.9700017,
        0.8500022,
        0.9699995,
        0.8499989,
        0.9586623,
        0.969997,
    ]

    Ptot = (2*GenP - LoadP) / base_power

    #LoadQ = LoadP .* sin.(acos.(cosϕ)) ./ cosϕ
    Qtot = Ptot .* sin.(acos.(cosϕ)) ./ cosϕ
    Qtot[2] = 0.0

    # Droop Control Parameters
    τ_Q = 8.0
    K_P = 10.0
    K_Q = 0.1
    V_r = 1.0
    τ_P =0.5

    begin
        buses = Array{PowerDynamics.AbstractNode,1}([])
        # push!(busses_static, SlackAlgebraic(U = 110E3 / base_voltage_HV))

        for (P, Q) in zip(Ptot, Qtot)
            push!(
                buses,
                VSIMinimal(τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q),
                #PQAlgebraic(S = - complex(P, Q) ./ base_power),
                )
        end
    end



    ###### Kanten ######



    # T = OLTC(
    #     from = 1,
    #     to = 2,
    #     Uos = 110E3 / base_voltage_HV, # Bemessungsspannung Oberspannungsseite in kV
    #     Uus = 20E3 / base_voltage,# Bemessungsspannung Unterspannungsseite in kV
    #     k = 0, # Kennzahl der Schaltgruppe, Bsp: YD5-> 5
    #     ssp = 10.0, # Stufenschalterposition
    #     stufe_us = 0.625, # Stufung pro Stufenschalterstellung in %
    #     Sr = 25E6 / base_power, #Bemessungsscheinleistung in MVA
    #     uk = 12.0, # Kurzschlussspannung in %
    #     Pvk = 25E3 / base_power, # Kupferverluste in MW
    #     Pvl = 0.0, # Eisenverluste in kW
    #     iLeer = 0.0, # Leerlaufstrom in %
    # )

    #ldata = CSV.read("$dir/lines.csv"; header = true)

    ldata = [ #L_km(nicht ändern)
        2.82,
        4.42,
        0.61,
        1.3,
        0.56,
        1.54,
        1.67,
        0.32,
        0.77,
        0.33,
        #0.24,
        #0.49,
    ]

    elist = [
        (1, 2),
        (2, 3),
        (3, 4),
        (3, 8),
        (4, 5),
        (5, 6),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 11),
        #(6, 7), # Schalter
        #(4, 11), # Schalter
    ]

    if close_switches == true
        append!(ldata,[0.24,0.49])
        append!(elist,[(6, 7),(4, 11)])
    end

    R1 = 0.501 # Ω/km
    X1 = 0.716 # Ω/km
    C1 = 0.1511749E-6 # F/km

    Zs = complex(R1, X1) .* ldata #[!, Symbol("L_km(nicht ändern)")]
    Yshs = 1im .* ω .* C1 .* ldata #[!, Symbol("L_km(nicht ändern)")]

    begin
        lines = Array{AbstractLine,1}([])
        # push!(lines, T)
        for (e, Z, Ysh) in zip(elist, Zs, Yshs)
            push!(
                lines,
                # PiModelLine(
                #     from = first(e),
                #     to = last(e),
                #     y = inv(Z) / base_admittance,
                #     y_shunt_km = Ysh / 2.0 / base_admittance,
                #     y_shunt_mk = Ysh / 2.0 / base_admittance,
                # ),
                StaticLine(
                    from = first(e),
                     to = last(e),
                     Y = inv(Z) / base_admittance,
                )
            )
        end
    end

    buses, lines, ldata
end