using Plots
const G = 9.80665

function get_dstate(s, μ)
    ṡ = copy(s)

    r = s[1]
    θ = s[2]
    vr = s[3]
    vθ = s[4]

    ṡ[1] = vr
    ṡ[2] = vθ
    ṡ[3] = (r*vθ^2.0 + G*(cos(θ) - μ)) / (1.00 + μ)
    ṡ[4] = - (2.00*vr*vθ + G*sin(θ)) / r
    return ṡ
end

function RK4!(motion, init_state, dt, μ)
    for i in 1:size(motion)[2]
        motion[:,i] = init_state

        ds1 = dt * get_dstate(init_state, μ)
        ds2 = dt * get_dstate(init_state + 0.5*ds1, μ)
        ds3 = dt * get_dstate(init_state + 0.5*ds2, μ)
        ds4 = dt * get_dstate(init_state + ds3, μ)
        init_state += (ds1 + 2*ds2 + 2*ds3 + ds4) / 6

        if init_state[1] <= 0.00 || init_state[1] >= 7.00
            break
        end
    end
end

function plotter(motion, time, μ)
    r = motion[1,:]
    θ = motion[2,:]
    x = r.*sin.(θ)
    y = - r.*cos.(θ)

    # plt1 = plot(time, motion[1,:]
    #     ,label="r-position"
    #     ,linecolor=:green
    # )
    
    # plot!(plt1, time, motion[3,:]
    #     ,label="r-velocity"
    #     ,linecolor=:maroon
    # )

    # plt2 = plot(time, motion[2,:]
    #     ,label="theta-position"
    #     ,linecolor=:green
    # )
    
    # plot!(plt2, time, motion[4,:]
    #     ,label="theta-velocity"
    #     ,linecolor=:maroon
    # )

    plt3 = plot(x, y
                ,label="μ = $(μ)"
                ,linecolor=:green
    )
    plot!(plt3
        ,xaxis="x-position"
        ,yaxis="y-position"
        ,legend=:bottomright
        ,aspect_ratio=:equal
    )
    
    # display(plt1);
    # display(plt2);
    display(plt3);
end

function SwingingAtwoodsMachine(μ::Float64=1.1185, t_fin::Float64=12.500)
    dt = 10.0^(-3.0)
    time = 0.0:dt:t_fin
    sys_motion = zeros(4, length(time))

    init_state = [1.00, pi/2.00, 0.00, 0.00]
    RK4!(sys_motion, init_state, dt, μ)

    plotter(sys_motion, time, μ)
end

function benchmark_func(μ::Float64=1.1185, t_fin::Float64=12.500)
    dt = 10.0^(-3.0)
    time = 0.0:dt:t_fin
    sys_motion = zeros(4, length(time))

    init_state = [1.00, pi/2.00, 0.00, 0.00]
    RK4!(sys_motion, init_state, dt, μ)
end