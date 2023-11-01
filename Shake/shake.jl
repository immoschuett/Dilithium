module SHAK3
export shake256
export sha3_256
export sha3_512
export shake128
# SHAKE 256 context structures
abstract type SHA end
abstract type SHAKE <: SHA end 
abstract type SHA3 <: SHA end 
mutable struct SHAKE_128_CTX <: SHAKE
    state::Array{UInt64,1}
    bytecount::UInt128
    buffer::Array{UInt8,1}
    bc::Array{UInt64,1}
    used::Bool
end
mutable struct SHAKE_256_CTX <: SHAKE
    state::Array{UInt64,1}
    bytecount::UInt128
    buffer::Array{UInt8,1}
    bc::Array{UInt64,1}
    used::Bool
end
mutable struct SHA3_512_CTX <: SHA3
    state::Array{UInt64,1}
    bytecount::UInt128
    buffer::Array{UInt8,1}
    bc::Array{UInt64,1}
    used::Bool
end
mutable struct SHA3_256_CTX <: SHA3
    state::Array{UInt64,1}
    bytecount::UInt128
    buffer::Array{UInt8,1}
    bc::Array{UInt64,1}
    used::Bool
end
# Round constants for SHA3 rounds
const AbstractBytes = Union{AbstractVector{UInt8},NTuple{N,UInt8} where N}
const SHA3_ROUND_CONSTS = UInt64[
    0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
    0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
    0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
    0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
    0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
    0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
    0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
    0x8000000000008080, 0x0000000080000001, 0x8000000080008008
]
# Rotation constants for SHA3 rounds
const SHA3_ROTC = UInt64[
    1,  3,  6,  10, 15, 21, 28, 36, 45, 55, 2,  14,
    27, 41, 56, 8,  25, 43, 62, 18, 39, 61, 20, 44
]
# Permutation indices for SHA3 rounds (+1'ed so as to work with julia's 1-based indexing)
const SHA3_PILN = Int[
    11, 8,  12, 18, 19, 4, 6,  17, 9,  22, 25, 5,
    16, 24, 20, 14, 13, 3, 21, 15, 23, 10,  7,  2
]

digestlen(::Type{SHAKE_128_CTX}) = 32
digestlen(::Type{SHA3_256_CTX}) = 32
digestlen(::Type{SHAKE_256_CTX}) = 32
digestlen(::Type{SHA3_512_CTX}) = 64
blocklen(::Type{SHAKE_128_CTX}) = UInt64(25*8 - 2*digestlen(SHAKE_128_CTX))
blocklen(::Type{SHA3_256_CTX}) = UInt64(25*8 - 2*digestlen(SHA3_256_CTX))
blocklen(::Type{SHAKE_256_CTX}) = UInt64(25*8 - 2*digestlen(SHAKE_256_CTX))
blocklen(::Type{SHA3_512_CTX}) = UInt64(25*8 - 2*digestlen(SHA3_512_CTX))
buffer_pointer(ctx::T) where {T<:SHAKE} = Ptr{state_type(T)}(pointer(ctx.buffer))
buffer_pointer(ctx::T) where {T<:SHA3} = Ptr{state_type(T)}(pointer(ctx.buffer))
# 64-bit Rotate-left (used in SHA3)
lrot(b,x,width) = ((x << b) | (x >> (width - b)))
L64(b,x) = lrot(b,x,64)
# construct an empty SHA context
SHAKE_128_CTX() = SHAKE_128_CTX(zeros(UInt64, 25), 0, zeros(UInt8, blocklen(SHAKE_128_CTX)), Vector{UInt64}(undef, 5), false)
SHAKE_256_CTX() = SHAKE_256_CTX(zeros(UInt64, 25), 0, zeros(UInt8, blocklen(SHAKE_256_CTX)), Vector{UInt64}(undef, 5), false)
SHA3_512_CTX()  = SHA3_512_CTX(zeros(UInt64, 25), 0, zeros(UInt8, blocklen(SHA3_512_CTX)), Vector{UInt64}(undef, 5), false)
SHA3_256_CTX()  = SHA3_256_CTX(zeros(UInt64, 25), 0, zeros(UInt8, blocklen(SHA3_256_CTX)), Vector{UInt64}(undef, 5), false)

function transform!(context::T) where {T<:SHA3}
    # First, update state with buffer
    pbuf = Ptr{eltype(context.state)}(pointer(context.buffer))
    for idx in 1:div(blocklen(T),8)
        context.state[idx] = context.state[idx] ⊻ unsafe_load(pbuf, idx)
    end
    bc = context.bc
    state = context.state
    # We always assume 24 rounds
    @inbounds for round in 0:23
        # Theta function
        for i in 1:5
            bc[i] = state[i] ⊻ state[i + 5] ⊻ state[i + 10] ⊻ state[i + 15] ⊻ state[i + 20]
        end
        for i in 0:4
            temp = bc[rem(i + 4, 5) + 1] ⊻ L64(1, bc[rem(i + 1, 5) + 1])
            j = 0
            while j <= 20
                state[Int(i + j + 1)] = state[i + j + 1] ⊻ temp
                j += 5
            end
        end
        # Rho Pi
        temp = state[2]
        for i in 1:24
            j = SHA3_PILN[i]
            bc[1] = state[j]
            state[j] = L64(SHA3_ROTC[i], temp)
            temp = bc[1]
        end
        # Chi
        j = 0
        while j <= 20
            for i in 1:5
                bc[i] = state[i + j]
            end
            for i in 0:4
                state[j + i + 1] = state[j + i + 1] ⊻ (~bc[rem(i + 1, 5) + 1] & bc[rem(i + 2, 5) + 1])
            end
            j += 5
        end
        # Iota
        state[1] = state[1] ⊻ SHA3_ROUND_CONSTS[round+1]
    end
    return context.state
end
function transform!(context::T) where {T<:SHAKE}
    # First, update state with buffer
    pbuf = Ptr{eltype(context.state)}(pointer(context.buffer))
    # after SHAKE_256_MAX_READ is reached, simply work with context.state[idx]
    if !context.used 
        for idx in 1:div(blocklen(T),8)
            context.state[idx] = context.state[idx] ⊻ unsafe_load(pbuf, idx)
        end
    end 
    bc = context.bc
    state = context.state
    # We always assume 24 rounds
    @inbounds for round in 0:23
        # Theta function
        for i in 1:5
            bc[i] = state[i] ⊻ state[i + 5] ⊻ state[i + 10] ⊻ state[i + 15] ⊻ state[i + 20]
        end
        for i in 0:4
            temp = bc[rem(i + 4, 5) + 1] ⊻ L64(1, bc[rem(i + 1, 5) + 1])
            j = 0
            while j <= 20
                state[Int(i + j + 1)] = state[i + j + 1] ⊻ temp
                j += 5
            end
        end
        # Rho Pi
        temp = state[2]
        for i in 1:24
            j = SHA3_PILN[i]
            bc[1] = state[j]
            state[j] = L64(SHA3_ROTC[i], temp)
            temp = bc[1]
        end
        # Chi
        j = 0
        while j <= 20
            for i in 1:5
                bc[i] = state[i + j]
            end
            for i in 0:4
                state[j + i + 1] = state[j + i + 1] ⊻ (~bc[rem(i + 1, 5) + 1] & bc[rem(i + 2, 5) + 1])
            end
            j += 5
        end
        # Iota
        state[1] = state[1] ⊻ SHA3_ROUND_CONSTS[round+1]
    end
    return context.state
end
function update!(context::T, data::U, datalen=length(data)) where {T<:SHA, U<:AbstractBytes}
    context.used && error("Cannot update CTX after `digest!` has been called on it")
    # We need to do all our arithmetic in the proper bitwidth
    UIntXXX = typeof(context.bytecount)
    # Process as many complete blocks as possible
    0 ≤ datalen ≤ length(data) || throw(BoundsError(data, firstindex(data)+datalen-1))
    len = convert(UIntXXX, datalen)
    data_idx = convert(UIntXXX, firstindex(data)-1)
    usedspace = context.bytecount % blocklen(T)
    while len - data_idx + usedspace >= blocklen(T)
        # Fill up as much of the buffer as we can with the data given us
        copyto!(context.buffer, usedspace + 1, data, data_idx + 1, blocklen(T) - usedspace)
        transform!(context)
        context.bytecount += blocklen(T) - usedspace
        data_idx += blocklen(T) - usedspace
        usedspace = convert(UIntXXX, 0)
    end
    # There is less than a complete block left, but we need to save the leftovers into context.buffer:
    if len > data_idx
        copyto!(context.buffer, usedspace + 1, data, data_idx + 1, len - data_idx)
        context.bytecount += len - data_idx
    end
end
# Finalize data in the buffer, append total bitlength, and return our precious hash!
function digest!(context::T) where {T<:SHA3}
    if !context.used
        usedspace = context.bytecount % blocklen(T)
        # If we have anything in the buffer still, pad and transform that data
        if usedspace < blocklen(T) - 1
            # Begin padding with a 0x06
            context.buffer[usedspace+1] = 0x06
            # Fill with zeros up until the last byte
            context.buffer[usedspace+2:end-1] .= 0x00
            # Finish it off with a 0x80
            context.buffer[end] = 0x80
        else
            # Otherwise, we have to add on a whole new buffer just for the zeros and 0x80
            context.buffer[end] = 0x06
            transform!(context)
            context.buffer[1:end-1] .= 0x0
            context.buffer[end] = 0x80
        end

        # Final transform:
        transform!(context)
        context.used = true
    end
    # Return the digest
    return reinterpret(UInt8, context.state)[1:digestlen(T)]
end
function digest!(context::T,d::UInt) where {T<:SHAKE}
    usedspace = context.bytecount % blocklen(T)
    # If we have anything in the buffer still, pad and transform that data
    if usedspace < blocklen(T) - 1
        # Begin padding with a 0x1f
        context.buffer[usedspace+1] = 0x1f
        # Fill with zeros up until the last byte
        context.buffer[usedspace+2:end-1] .= 0x00
        # Finish it off with a 0x80
        context.buffer[end] = 0x80
    else
        # Otherwise, we have to add on a whole new buffer just for the zeros and 0x80
        context.buffer[end] = 0x1f
        transform!(context)
        context.buffer[1:end-1] .= 0x0
        context.buffer[end] = 0x80
    end
    # Final transform:
    transform!(context)
    # Return the digest
    if d <= blocklen(T)

        return reinterpret(UInt8, context.state)[1:d]
    else 
        a = reinterpret(UInt8, context.state)[1:blocklen(T)]
        context.used = true
        b = digest!(context,d-blocklen(T))
        return vcat(a,b)
    end
end
function sha3_256(data::AbstractBytes)
    ctx = SHA3_256_CTX()
    update!(ctx, data)
    return digest!(ctx)
end
function sha3_512(data::AbstractBytes)
    ctx = SHA3_512_CTX()
    update!(ctx, data)
    return digest!(ctx)
end
function shake128(data::AbstractBytes,d::UInt)
    ctx = SHAKE_128_CTX()
    update!(ctx, data)
    return digest!(ctx,d)
end
function shake256(data::AbstractBytes,d::UInt)
    ctx = SHAKE_256_CTX()
    update!(ctx, data)
    return digest!(ctx,d)
end
end # module
