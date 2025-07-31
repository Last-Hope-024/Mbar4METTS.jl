using FFTW
using Statistics

function get_equilibration(data::AbstractVector{Float64})
    """ Detect equilibrated region of a dataset by maximizing the number of effectively uncorrelated samples.
    """
    N = length(data)
    if std(data) == 0.0
        return 1, 1, N
    end
    gt = ones(Float64,N-1)
    Neff_t = ones(Float64, N-1)
    
for t in 1:(N-1)
    gt[t] = compute_integrated_autocorr(data[t:end])
    Neff_t[t] = (N - t) / gt[t]  # effective sample size
end
    
    Neff_max = maximum(Neff_t)
    t_max = argmax(Neff_t)
    g = gt[t_max]
    t_max, g, Neff_max
    
    array = collect(range(t_max, step=g, stop=N))
    array = round.(Int,array)
    return array
end


function compute_integrated_autocorr(data::AbstractVector{Float64})
    """
    Compute integrated autocorrelation time
    """
    N = length(data)
    Ct = compute_auto_corr(data)[1:end]
    t_grid = collect(1:N)
    g_t = 2.0 .*Ct .*(1.0 .- t_grid/N)
    filter!(>(0.0),g_t)
    return max(1.0,1.0 + sum(g_t[2:end])) # carefull avoid negative values
end

function compute_auto_corr(data::AbstractVector{Float64})
    """
    Compute autocorrelation time
    """
    N = length(data)
    data2 = data .- mean(data)
    n_padded = 2^ceil(Int, log2(2*N - 1))  # minimal safe padding
    fx = fft([data2; zeros(n_padded - N)])
    acf = real(ifft(abs2.(fx)))[1:N]
    acf ./= acf[1]  # Normalise by variance
    return acf
end