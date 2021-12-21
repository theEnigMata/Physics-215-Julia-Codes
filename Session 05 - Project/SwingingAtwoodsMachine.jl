using Plots
const G = 9.80665

function get_dstate(s, mu)
    ds = copy(s)

    r = s[1]
    theta = s[2]
    vr = s[3]
    vtheta = s[4]

    ds[1] = vr
    ds[2] = vtheta
    ds[3] = (r*vtheta^2 + G*(cos(theta) - mu)) / (1.00 + mu)
    ds[4] = - (2.00*vr*vtheta + G*sin(theta)) / r
    return ds
end

function RK4!(motion, init_state, dt, mu)
    motion[:,1] = init_state
    for i in 2:size(motion)[2]
        ds1 = dt * get_dstate(init_state, mu)
        ds2 = dt * get_dstate(init_state + 0.5*ds1, mu)
        ds3 = dt * get_dstate(init_state + 0.5*ds2, mu)
        ds4 = dt * get_dstate(init_state + ds3, mu)
        init_state += (ds1 + 2*ds2 + 2*ds3 + ds4) / 6
        motion[:,i] = init_state
        if init_state[1] <= 0.00
            break
        end
    end
end

function plot_them(motion, time)
    x = zeros(Float64, size(time))
    y = zeros(Float64, size(time))
    x .= motion[1,:].*sin.(motion[2,:])
    y .= - motion[1,:].*cos.(motion[2,:])

    plt1 = plot(time, motion[1,:]
        ,label="r-position"
        ,linecolor=:green
    )
    
    plot!(plt1, time, motion[2,:]
        ,label="r-velocity"
        ,linecolor=:maroon
    )


    plt3 = plot(x, y
                ,label="m-trajectory"
                ,linecolor=:green
    )
    plot!(plt3
        ,xaxis="x-position"
        ,yaxis="y-position"
        ,legend=:bottomright
        ,aspect_ratio=:equal
    )
    
    display(plt3);
end

function SwingingAtwoodsMachine(mu::Float64=1.1185, t_fin::Float64=12.500)
    dt = 10.0^(-2.00)
    time = 0:dt:t_fin
    sys_motion = zeros(4, length(time))

    init_state = [1.00, 3.0*pi/4.00, 0.00, 0.00]
    RK4!(sys_motion, init_state, dt, mu)

    plot_them(sys_motion, time)
end

