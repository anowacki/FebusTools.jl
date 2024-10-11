module FebusTools

using Dates: Dates, DateTime, Second, Millisecond, Microsecond, datetime2unix, unix2datetime
using HDF5: attrs, h5open
using Statistics: mean

export
    read_hdf5


"Element type of Febus A1 data"
const _A1_ELTYPE = Float32

struct FebusData{M,Da,Ti,Di}
    "Matrix of strain or strain rate with dimensions (time, distance)"
    data::M
    "Dates of each sample"
    dates::Da
    "Times of each sample"
    times::Ti
    "Distances of each channel"
    distances::Di
    "Other information read from file"
    metadata::Dict{Symbol,Any}
end

# Keyword argument constructor
FebusData(; data, dates, times, distances, metadata) =
    FebusData(data, dates, times, distances, metadata)

"""
    read_hdf5(file; xlim=(-Inf, Inf), xdecimate=1, tlim=($(typemin(DateTime)), $(typemax(DateTime))), blocks, zones=(==(1)), sources=(==(1)), version=nothing) -> ::FebusData

Read strain or strain rate data from `file` on disk.  By default, all data are
read into memory.

If `xlim` is given a length-2 collection, then only channels whose distances
lie between `xlim[begin]` and `xlim[end]` will be read.  Likewise, only blocks with
times between `tlim[begin]` and `tlim[end]` are read.  Note that whole blocks
(usually of duration 1 s) are always read, so the stand and end times of the
returned data can be respectively before and after `xlim[begin]` and `xlim[end]`.
To read only the `n`th channel, set `xdecimate=n`.

Instead of `tlim`, you can also choose which blocks are read in with the
`blocks` keyword argument, which again takes two items, the start and end
block index (which begins at 1).  Block indices outside the range of the data
are warned about and only existing data which match the range are read.

To limit data to a set of zones, supply a function of form `f(::Int)::Bool`
to `zones` which returns `true` when its first argument is the number of the
zone to accept.

Likewise, supply a function of the form `f(::Int)::Bool` for `sources` to accept
source a set of source numbers.

The data type (i.e., strain or strain rate) is recorded in the `:data_type`
field of the metada dictionary.

By default, the file version is determined from the HDF5 file attributes,
noting that files without a version attribute are considered to be 'v1'
files.  Pass a version (e.g., `v"2.1"`) to force the assumption of a
specific file version when parsing the data.

!!! note
    Currently only the first source of the first zone is read, and an error
    is thrown if more than one zone or source is present in a file or if
    anything other than the first zone is requested.  This will change in
    a breaking update to FebusTools.
"""
function read_hdf5(file;
    tlim::Union{Tuple,AbstractArray,Nothing}=nothing,
    xlim::Union{Tuple,AbstractArray}=(-Inf, Inf),
    xdecimate::Integer=1,
    blocks::Union{Tuple{Integer,Integer},AbstractArray{<:Integer},Nothing}=nothing,
    zones=(==(1)),
    sources=(==(1)),
    header_only::Bool=false,
    version::Union{Nothing,VersionNumber}=nothing,
)
    zones == (==(1)) && sources == (==(1)) ||
        throw(ArgumentError("zones and sources are not currently supported"))

    # Requirement to pass two items and not a range or filtering function is
    # because `HDF5` only slices with continuous ranges, which we can construct
    # from length-2 collections.
    if !isnothing(tlim) && length(tlim) != 2
        throw(ArgumentError("`tlim` must contain two items"))
    end
    if length(xlim) != 2
        throw(ArgumentError("`xlim` must contain two items"))
    end
    if xdecimate < 1
        throw(ArgumentError("`xdecimate` must be 1 or more"))
    end

    # Allow either time limiting by block indices or date, but not both
    if !isnothing(tlim) && !isnothing(blocks)
        throw(ArgumentError("only one of `tlim` or `blocks` can be used"))
    elseif isnothing(tlim) && isnothing(blocks)
        tlim = (typemin(DateTime), typemax(DateTime))
    end

    metadata = Dict{Symbol,Any}()

    h5open(file, "r") do f
        # Get sensor group
        sensor_key = only(keys(f["/"]))
        sensor = f["/"][sensor_key]
        # Get source
        source_key = only(keys(sensor))
        source = sensor[source_key]

        # Get file version, or use the assumed version passed in
        # There are at least two file versions in the wild; one without a
        # `"Version"` attribute which we call v1, and one with a
        # version attribute.  We have seen at least v2.3.21.
        file_version = if isnothing(version)
            v = VersionNumber(get(attrs(source), "Version", "1.0.0"))
            if v < v"1" || v >= v"3"
                @warn("File version number ($v) is outside expected range for file $file")
            end
            v
        else
            version
        end

        # Get zones
        zone_key = only(filter(!=("time"), keys(source)))
        zone = source[zone_key]

        # Global attributes
        # Total size of blocks before downsampling
        whole_extent = metadata[:WholeExtent] = Int.(attrs(source)["WholeExtent"])

        # Derived attributes using raw values
        zone_atts = attrs(zone)
        # Block overlap in percent
        block_overlap_percent = if file_version < v"2.3.13"
            # Original files always have complete overlap, per Table 2
            # of user guide v5.
            100
        elseif file_version >= v"2.3.13"
            # User guide v10 specifies variable overlap, to a minimum of
            # 5%
            Int(only(zone_atts["BlockOverlap"]))
        end
        metadata[:BlockOverlap] = block_overlap_percent
        # Block rate mHz and in-block sample spacing (ms)
        sampling_interval = zone_atts["Spacing"][2]/1000
        sampling_rate = 1000/zone_atts["Spacing"][2]

        # Get attributes for the first zone and convert to default types
        for (k, v) in zone_atts
            key = Symbol(k)
            val = length(v) == 1 ? only(v) : v

            # Spatial sampling resolution
            if key == :SamplingRes
                # Convert from cm -> m
                val /= 100
            # Rate of block writing
            elseif key == :BlockRate
                # Block rate is in mHz
                metadata[:BlockLength] = (100 + block_overlap_percent)*10/val
                metadata[:BlockInterval] = 1000/val
                # mHz -> Hz
                val /= 1000
            # Frequency of laser pulses
            elseif key == :PulseRateFreq
                # mHz -> Hz
                val = val/1000
            elseif key == :DataDomain
                val = Int.(val)
            # Extent of blocks for this zone
            elseif key == :Extent
                val = Int.(val)
            # Raw sampling rate of instrument
            elseif key == :SamplingRate
                # mHz -> Hz
                val /= 1_000
            elseif key == :DerivationTime
                # ms -> s
                val /= 1000
            end

            metadata[key] = val
        end

        # Timestamps of each block
        # For older files, the timestamp was the start of the block.
        # Files with version >= v2.3.13 have the middle of the block
        # as the time (per V10 of the manual)
        block_start_times = source["time"][:] .+ metadata[:Origin][2]/1000
        block_start_dates = unix2datetime.(block_start_times)

        # Number of points per whole block
        metadata[:nsamples_per_block] = whole_extent[4] - whole_extent[3] + 1
        if metadata[:nsamples_per_block] != metadata[:BlockLength]*metadata[:PulseRateFreq]
            @warn "Unexpected total block extent ($whole_extent) for block " *
                "length ($(metadata[:BlockLength]) s) and pulse rate frequency " *
                "($(metadata[:PulseRateFreq]) Hz)"
        end

        # Find distances along the fibre
        x0::Float64 = metadata[:Origin][1]
        distances = x0 .+ (metadata[:Extent][1]:metadata[:Extent][2]).*metadata[:Spacing][1]

        # Data.  Take strain if present, otherwise strain rate.
        raw_data = if file_version < v"2.3.13"
            if haskey(zone, "Strain")
                metadata[:data_type] = "strain [nstrain]"
                zone["Strain"]
            elseif haskey(zone, "StrainRate")
                metadata[:data_type] = "strain rate [nstrain/s]"
                zone["StrainRate"]
            else
                error("no known data fields present in file $file")
            end
        else
            if haskey(zone, "Strain [nStrain]")
                metadata[:data_type] = "strain [nstrain]"
                zone["Strain [nStrain]"]
            elseif haskey(zone, "Strain Rate [nStrain|s]")
                metadata[:data_type] = "strain rate [nstrain/s]"
                zone["Strain Rate [nStrain|s]"]
            else
                error("no known data fields present in file $file")
            end
        end

        # Calculate which blocks to include
        total_nblocks = size(raw_data, 3)

        # Just return metadata if requested
        if header_only
            return FebusData(Matrix{_A1_ELTYPE}(undef, 0, 0), [], [], distances, metadata)
        end

        if !isnothing(tlim)
            cut1, cut2 = tlim
            t_ind1 = findfirst(x -> cut1 <= x <= cut2, block_start_dates)
            t_ind2 = findlast(x -> cut1 <= x <= cut2, block_start_dates)
            if isnothing(t_ind1)
                @warn "No blocks between $cut1 and $cut2 in file '$file'"
                return FebusData(Matrix{_A1_ELTYPE}(undef, 0, 0), [], [], distances, metadata)
            end
        else
            t_ind1, t_ind2 = blocks
            if t_ind1 < 1
                @warn "Start block $(t_ind1) is below data range (i.e., 1) in file '$file'"
            end
            if t_ind2 > total_nblocks
                @warn "End block $(t_ind2) is after the last block ($(total_nblocks)) in file '$file'"
            end
        end

        nblocks = t_ind2 - t_ind1 + 1

        # Calculate which channels to include
        x1, x2 = xlim
        x_ind1 = findfirst(x -> x1 <= x <= x2, distances)
        x_ind2 = findlast(x -> x1 <= x <= x2, distances)
        if isnothing(x_ind1)
            @warn "No channels between $x1 and $x2 m in file '$file'"
            return FebusData(Matrix{_A1_ELTYPE}(undef, 0, 0), [], [], distances, metadata)
        end
        nchannels = length(x_ind1:xdecimate:x_ind2)

        # Calculate which block sample numbers we need to avoid overlap (starting at 0)
        block_num1, block_num2 = metadata[:Extent][3], metadata[:Extent][4]
        ninds = block_num2 - block_num1 + 1
        if ninds*sampling_interval < metadata[:BlockInterval]
            error("block extent ($(block_num1:block_num2)) does not cover whole ",
                "block spacing ($(metadata[:BlockInterval]) s).  Cannot read ",
                "continuous data from file '$file'.")
        end
        block_num2 = block_num1 + round(Int, metadata[:BlockInterval]/sampling_interval) - 1

        # Convert to Julia indices of HDF5 array block samples (starting at 1),
        # always starting with the first available sample.
        # TODO: Check if a later sample aligns better with UTC and use that?
        block_ind1 = 1
        block_ind2 = block_num2 - block_num1 + 1

        nsamples = length(block_num1:block_num2)*nblocks

        # Cut out data
        # Reading becomes slow when xdecimate is small, because
        # it seems HDF5 makes lots of reads (perhaps one for each
        # slice of the first dimension).  If `xdecimate`` is small,
        # then chunk the reads up along a more sensible dimension,
        # and decimate the first afterwards.
        if xdecimate == 1 || xdecimate >= 50
            unshaped_data =
                raw_data[x_ind1:xdecimate:x_ind2, block_ind1:block_ind2, t_ind1:t_ind2]
            data = permutedims(reshape(unshaped_data, nchannels, :))
        else
            # Avoid slowdown for small decimation values
            # by just reading it all in then removing the channels
            # we don't want.
            # FIXME: Implement a better time-chunked strategy as this
            #        requires us to read it all into memory.
            unshaped_data =
                raw_data[x_ind1:x_ind2, block_ind1:block_ind2, t_ind1:t_ind2]
            data = permutedims(reshape(unshaped_data[1:xdecimate:end,:,:],
                nchannels, :))
        end

        if size(data) != (nsamples, nchannels)
            error("error in calculating data slice.  Expected shape ",
                (nsamples, nchannels), "; got ", size(data))
        end

        block_dates = block_start_dates[t_ind1:t_ind2]
        times = (block_num1 .+ (0:(nsamples - 1))) .* sampling_interval
        cut_distances = distances[x_ind1:xdecimate:x_ind2]

        FebusData(data, block_dates, times, cut_distances, metadata)
    end
end

end # module
