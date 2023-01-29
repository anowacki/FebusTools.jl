module FebusTools

using Dates: Dates, DateTime, Second, Millisecond, datetime2unix, unix2datetime
using HDF5: attrs, h5open
using NanoDates: NanoDate
using Statistics: mean

export
    read_hdf5


struct FebusData{M,Da,Ti,Di}
    data::M
    dates::Da
    times::Ti
    distances::Di
    metadata::Dict{Symbol,Any}
end

"""
    read_hdf5(file[, cut1[, cut2]]; zones, sources) -> ::FebusData

Read strain or strain rate data from `file` on disk, reading blocks between
dates `cut1` and `cut2`.  By default, all data are read.

To limit data to a set of zones, supply a function of form `f(::Int)::Bool`
to `zones` which returns `true` when its first argument is the number of the
zone to accept.

Likewise, supply a function of the form `f(::Int)::Bool` for `sources` to accept
source a set of source numbers.

!!! note
    Currently only the first source of the first zone is read, and an error
    is thrown if more than one zone or source is present in a file.
"""
function read_hdf5(file, cut1=typemin(DateTime), cut2=typemax(DateTime);
    zones=(>=(1)), sources=(>=(1))
)
    metadata = Dict{Symbol,Any}()

    h5open(file, "r") do f
        # Get sensor group
        sensor_key = only(keys(f["/"]))
        sensor = f["/"][sensor_key]
        # Get source
        source_key = only(keys(sensor))
        source = sensor[source_key]
        # Get times
        block_times = source["time"][:]
        block_dates = unix2datetime.(block_times)
        # Get zones
        zone_key = only(filter(!=("time"), keys(source)))
        zone = source[zone_key]

        # Get attributes for the first zone and convert to default types
        for (k, v) in attrs(zone)
            key = Symbol(k)
            val = length(v) == 1 ? only(v) : v

            if key == :SamplingRes
                # Convert from cm -> m
                val /= 100
            elseif key == :BlockRate
                metadata[:BlockLength] = 2000/val
                # mHz -> Hz
                val /= 1000
            elseif key == :PulseRateFreq
                # mHz -> Hz
                val = Int(val)/1000
            elseif key == :DataDomain
                val = Int.(val)
            elseif key == :Extent
                val = Int.(val)
            elseif key == :SamplingRate
                metadata[:SamplingInterval] = Dates.Microsecond(1_000_000_000_000//val)
                # µHz -> Hz
                val /= 1_000_000
            elseif key == :DerivationTime
                # ms -> s
                val /= 1000
            end

            metadata[key] = val
        end

        # Calculate recorded sampling interval in s
        sampling_interval = metadata[:BlockLength]/metadata[:Extent][2]

        # Number of points per whole block
        metadata[:nsamples_per_block] = metadata[:BlockLength]*metadata[:SamplingRate]

        x0::Float64 = metadata[:Origin][1]
        distances = x0 .+ (metadata[:Extent][1]:metadata[:Extent][2]).*metadata[:Spacing][1]
        # distances = 0:metadata[:Spacing][1]:metadata[:FiberLength]

        ind1 = findfirst(x -> cut1 <= x <= cut2, block_dates)
        ind2 = findlast(x -> cut1 <= x <= cut2, block_dates)
        if isnothing(ind1)
            @warn "No blocks between $cut1 and $cut2 in file '$file'"
            return FebusData(Matrix{Float32}(undef, 0, 0), [], [], distances, metadata)
        end
        nblocks = ind2 - ind1 + 1
        data = let d = zone["StrainRate"]
            permutedims(reshape(d[:,:,ind1:ind2], size(d, 1), :))
        end
        dates = block_dates[ind1] .+
            (Millisecond(0):metadata[:SamplingInterval]:Millisecond(1000*nblocks))
        times = (0:(size(data, 1) - 1)) .* sampling_interval

        FebusData(data, dates, times, distances, metadata)
    end
end

end # module
