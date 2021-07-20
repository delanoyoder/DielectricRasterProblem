module DielectricRaster

using StaticArrays

struct Square{T}
    x::Float64
    y::Float64
    w::Float64
    dielectric::T
end

struct GriddedArray{T}
    xlims::Tuple{Float64, Float64}
    ylims::Tuple{Float64, Float64}
    pixels::Matrix{T}
end

const TensorGriddedArray = GriddedArray{SMatrix{2, 2, Float64}}

## Q1
function raster_area!(A::GriddedArray, s::Square)

    # Locating the lower left and upper right vertices of the Square.
    square_ll = (s.x-s.w/2, s.y-s.w/2)
    square_ur = (s.x+s.w/2, s.y+s.w/2)

    # Storing the width, height, and area of a pixel to a variable for readability.
    pixel_w = A.xlims[2] / (size(A.pixels)[1] - 1)
    pixel_h = A.ylims[2] / (size(A.pixels)[2] - 1)
    pixel_a = pixel_w * pixel_h

    # Looping through the pixel array by the pixel centers.
    for x=A.xlims[1] : pixel_w : A.xlims[2], 
        y=A.ylims[1] : pixel_h : A.ylims[2]

        # Locating the lower left and upper right vertices of the pixel.
        pixel_ll = (x - pixel_w/2, y - pixel_h/2)
        pixel_ur = (x + pixel_w/2, y + pixel_h/2)
        
        # Calculating the overlapping area of the pixel and Square.
        ol_w = min(pixel_ur[1], square_ur[1]) - max(pixel_ll[1], square_ll[1])
        ol_h = min(pixel_ur[2], square_ur[2]) - max(pixel_ll[2], square_ll[2])

        # Calculating the current pixel indices.
        pixel_x = Int(x / pixel_w) + 1
        pixel_y = Int(y / pixel_h) + 1

        # If the pixel width or height is negative, the area is not overlapping and should be 0.
        if ol_w > 0 && ol_h > 0

            # Calculating the fractional overlap of the Square and pixel and saving it to the respective pixel value.
            ol_a = ol_w * ol_h
            frac_ol =  ol_a / pixel_a
            A.pixels[pixel_x, pixel_y] = A.pixels[pixel_x, pixel_y] * frac_ol
        else
            A.pixels[pixel_x, pixel_y] = 0.0
        end
    end

end

## Q2
function raster_scalar_dielectric!(A::GriddedArray, s::Square)

    # Determining the weights for the Square.
    weights = GriddedArray(A.xlims, A.ylims, ones(size(A.pixels)))
    raster_area!(weights, s)
    
    # Saving the weighted sums to the respective pixels in the GriddedArray object.
    A.pixels[:,:] = (weights.pixels * s.dielectric) + ((1.0 .- weights.pixels) .* A.pixels)
end

## Q3
function raster_harmonic!(A::GriddedArray, s::Square)

    # Determining the weights for the Square.
    weights = GriddedArray(A.xlims, A.ylims, ones(size(A.pixels)))
    raster_area!(weights, s)

    # Saving the harmonic average to the respective pixels in the GriddedArray object.
    A.pixels[:,:] = ones(size(A.pixels)) ./ ((weights.pixels ./ s.dielectric) + ((1.0 .- weights.pixels) ./ A.pixels))

end

## Q4
function find_nearest_surface_normal(s::Square, coord)

    # This function uses the assumption that points on a diagonal of the Square are ill-defined.
    # If a point lies on a diagonal, the nearest vertex is used as the closest element of the Square.

    # Translating the Square's center to (0,0) and combating rounding errors.
    # (Not sure what the best practice is for correcting rounding errors in Julia)
    x,y = round(coord[1] - s.x, digits=15), 
          round(coord[2] - s.y, digits=15)

    # Checking for points on the line y=x.
    if y == x
        if x > 0
            return SVector(1/sqrt(2), 1/sqrt(2))
        else
            return SVector(-1/sqrt(2), -1/sqrt(2))
        end
    end

    # Checking for points on the line y=-x.
    if y == -x
        if x > 0
            return SVector(1/sqrt(2), -1/sqrt(2))
        else
            return SVector(-1/sqrt(2), 1/sqrt(2))
        end
    end

    # Checking points nearest an edge.
    if x > -s.w/2 && x < s.w/2 || y > -s.w/2 && y < s.w/2

        # Checking points outside the Square.
        if x > s.w/2
            return SVector(1.0, 0.0)
        elseif y > s.w/2
            return SVector(0.0, 1.0)
        elseif x < -s.w/2
            return SVector(-1.0, 0.0)
        elseif y < -s.w/2
            return SVector(0.0, -1.0)
        end

        # Checking points inside the Square.
        d = [abs(s.w/2 - x),
             abs(s.w/2 - y),
             abs(-s.w/2 - x),
             abs(-s.w/2 - y)]
        vects = [SVector(1.0, 0.0),
                 SVector(0.0, 1.0),
                 SVector(-1.0, 0.0),
                 SVector(0.0, -1.0)]
        return vects[argmin(d)]
    end

    # Checking points nearest a vertex.
    d = [sqrt(.5 * s.w^2 - s.w * (x + y) + x^2 + y^2),
         sqrt(.5 * s.w^2 + s.w * (x - y) + x^2 + y^2),
         sqrt(.5 * s.w^2 - s.w * (x - y) + x^2 + y^2),
         sqrt(.5 * s.w^2 + s.w * (x + y) + x^2 + y^2)]
    vects = [SVector(1/sqrt(2), 1/sqrt(2)),
             SVector(-1/sqrt(2), 1/sqrt(2)),
             SVector(1/sqrt(2), -1/sqrt(2)),
             SVector(-1/sqrt(2), -1/sqrt(2))]
    return vects[argmin(d)]

end

## Q5
function orient_along_normal(A, normal)

    # Calculating the inverse of the new basis.
    basis = inv([normal[2] normal[1]; -normal[1] normal[2]])

    # Applying the inverse basis to our vector or matrix.
    return basis * A
end

## Q6, OPTIONAL
function raster_tensor_dielectric!(A::GriddedArray{<:SMatrix}, B::GriddedArray, s::Square)
    # implement me
end

export Square, GriddedArray, TensorGriddedArray
export raster_area!, raster_scalar_dielectric!, raster_harmonic!
export find_nearest_surface_normal, orient_along_normal, raster_tensor_dielectric!

end
