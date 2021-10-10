function mandelbrot_iterations(c::Complex, max_iterations::Int64)
    z = complex(0.0, 0.0)
    for iteration=1:max_iterations
        z = z^2 + c
        if abs2(z) >= 4
            return iteration
        end
    end
    return 0
end

function mandelbrot(a::Float64=0.00, b::Float64=0.00, n::Int64=3500)
    if a == 0 && b == 0
        re_min, re_max = -1.50, 0.50
        im_min, im_max = -1.25, 1.25
    else
        re_min, re_max = a - 0.125, a + 0.125
        im_min, im_max = b - 0.125, b + 0.125
    end
    
    complex_set    = zeros(ComplexF64, (n, n))
    mandelbrot_set = zeros(Int64, (n, n))
    
    y = 1
    for z_re = range(re_min, re_max, length=n)
        x = 1
        for z_im = range(im_min, im_max, length=n)
            complex_set[x, y] = complex(z_re, z_im)
            x += 1
        end
        y += 1
    end
        
    mandelbrot_set .= mandelbrot_iterations.(complex_set, 250)
    
    heatmap(mandelbrot_set, c=:inferno, axis=nothing, aspect_ratio=:equal, legend=:none, dpi=480, fmt=:png)
end;