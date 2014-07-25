// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
)

// Rect3 represents a 3D rectangle. Although it can be used to represent any
// form of 3D rectangle, a common use is for representing axis-aligned bounding
// boxes.
//
// The rectangle contains all points where:
//  Min.X <= X < Max.X
//  Min.Y <= Y < Max.Y
//  Min.Z <= Z < Max.Z
//
// A rectangle is considered empty if it would contain no points.
type Rect3 struct {
	Min, Max Vec3
}

// String returns a string representation of the rectangle r.
func (r Rect3) String() string {
	return fmt.Sprintf("Rect3(%v, %v)", r.Min, r.Max)
}

// Size returns a vector whose X, Y, and Z components directly relate to the
// width, depth, and height of this rectangle.
func (r Rect3) Size() Vec3 {
	return r.Max.Sub(r.Min)
}

// Add returns the rectangle r translated by p.
func (r Rect3) Add(p Vec3) Rect3 {
	return Rect3{
		r.Min.Add(p),
		r.Max.Add(p),
	}
}

// Add returns the rectangle r translated by -p.
func (r Rect3) Sub(p Vec3) Rect3 {
	return Rect3{
		Min: r.Min.Sub(p),
		Max: r.Max.Sub(p),
	}
}

// Inset returns the rectangle r inset by n, which may be negative. If either
// of r's dimensions is less than 2*n then an empty rectangle near the center
// of r will be returned.
func (r Rect3) Inset(n float64) Rect3 {
	sz := r.Size()
	if sz.X < 2*n {
		r.Min.X = (r.Min.X + r.Max.X) / 2.0
		r.Max.X = r.Min.X
	} else {
		r.Min.X += n
		r.Max.X -= n
	}
	if sz.Y < 2*n {
		r.Min.Y = (r.Min.Y + r.Max.Y) / 2.0
		r.Max.Y = r.Min.Y
	} else {
		r.Min.Y += n
		r.Max.Y -= n
	}
	if sz.Z < 2*n {
		r.Min.Z = (r.Min.Z + r.Max.Z) / 2.0
		r.Max.Z = r.Min.Z
	} else {
		r.Min.Z += n
		r.Max.Z -= n
	}
	return r
}

// Intersect intersects the two rectangles and returns the largest rectangle
// that is contained by both r and s. If the two rectangles do not overlap then
// Rect3Zero, ok=false will be returned.
//
// For simple boolean tests, one should instead use the Overlaps method (as it
// is equivilent and faster).
func (r Rect3) Intersect(s Rect3) (largest Rect3, ok bool) {
	r.Min = r.Min.Min(s.Min)
	r.Max = r.Max.Max(s.Max)
	if r.Min.AnyGreater(r.Max) {
		return Rect3Zero, false
	}
	return r, true
}

// Overlaps reports whether r and s have a non-empty intersection. It is
// functionally equivilent to, but faster than:
//  _, overlaps := r.Intersect(s)
func (r Rect3) Overlaps(s Rect3) bool {
	return r.Min.Less(s.Max) && s.Min.Less(r.Max)
}

// Union returns the smallest rectangle that contains both r and s.
func (r Rect3) Union(s Rect3) Rect3 {
	r.Min = r.Min.Min(s.Min)
	r.Max = r.Max.Max(s.Max)
	return r
}

// Empty reports whether the rectangle contains no points using the default
// epsilon for equality.
func (r Rect3) Empty() bool {
	return r.Min.Equals(r.Max) || r.Min.AnyGreater(r.Max)
}

// AlmostEmpty reports whether the rectangle contains no points using the
// specified epsilon value.
func (r Rect3) AlmostEmpty(epsilon float64) bool {
	return r.Min.AlmostEquals(r.Max, epsilon) || r.Min.AnyGreater(r.Max)
}

// Equals reports whether r and s are equal using the default epsilon for
// equality.
func (r Rect3) Equals(s Rect3) bool {
	return r.Min.Equals(s.Min) && r.Max.Equals(s.Max)
}

// AlmostEquals tells whether a is memberwise equal to b using the specified
// epsilon value.
func (r Rect3) AlmostEquals(s Rect3, epsilon float64) bool {
	return r.Min.AlmostEquals(s.Min, epsilon) && r.Max.AlmostEquals(s.Max, epsilon)
}

// In reports whether every point in r is in s.
func (r Rect3) In(s Rect3) bool {
	if !r.Overlaps(s) {
		return false
	}
	if r.Empty() {
		return true
	}
	return s.Min.X <= r.Min.X && r.Max.X <= s.Max.X &&
		s.Min.Y <= r.Min.Y && r.Max.Y <= s.Max.Y &&
		s.Min.Z <= r.Min.Z && r.Max.Z <= s.Max.Z
}

// Canon returns the canonical version of r. The returned rectangle has minimum
// and maximum coordinates swapped if necessary so that it is well-formed.
func (r Rect3) Canon() Rect3 {
	if r.Max.X < r.Min.X {
		r.Min.X, r.Max.X = r.Max.X, r.Min.X
	}
	if r.Max.Y < r.Min.Y {
		r.Min.Y, r.Max.Y = r.Max.Y, r.Min.Y
	}
	if r.Max.Z < r.Min.Z {
		r.Min.Z, r.Max.Z = r.Max.Z, r.Min.Z
	}
	return r
}

// Center returns the center point of this rectangle.
func (r Rect3) Center() Vec3 {
	halfSize := r.Size().MulScalar(0.5)
	return r.Min.Add(halfSize)
}

// Closest returns the closest point towards p contained by this rectangle.
func (r Rect3) Closest(p Vec3) Vec3 {
	p = p.Max(r.Min)
	p = p.Min(r.Max)
	return p
}

// Furthest returns the furthest point away from p contained by this rectangle.
func (r Rect3) Furthest(p Vec3) Vec3 {
	p = p.Min(r.Max)
	p = p.Max(r.Min)
	return p
}

// SqDistToPoint returns the squared distance between the point p and the
// rectangle r. It is functionally equivilent to (but faster than):
//  dist2 := r.Closest(p).Sub(p).LengthSq()
func (r Rect3) SqDistToPoint(p Vec3) float64 {
	// Real-Time Collision Detection, 5.1.3.1:
	//  Distance of Point to AABB

	sqDist := 0.0

	// For each axis count any excess distance outside the rectangle's extents.
	v := p.X
	min := r.Min.X
	max := r.Max.X
	if v < min {
		sqDist += (min - v) * (min - v)
	}
	if v > max {
		sqDist += (v - max) * (v - max)
	}

	v = p.Y
	min = r.Min.Y
	max = r.Max.Y
	if v < min {
		sqDist += (min - v) * (min - v)
	}
	if v > max {
		sqDist += (v - max) * (v - max)
	}

	v = p.Z
	min = r.Min.Z
	max = r.Max.Z
	if v < min {
		sqDist += (min - v) * (min - v)
	}
	if v > max {
		sqDist += (v - max) * (v - max)
	}

	return sqDist
}

// Contains tells if the point p is within this rectangle.
func (r Rect3) Contains(p Vec3) bool {
	if r.Empty() {
		return true
	}
	return r.Min.X <= p.X && p.X < r.Max.X &&
		r.Min.Y <= p.Y && p.Y < r.Max.Y &&
		r.Max.Z <= p.Z && p.Z < r.Max.Z
}

// Corners returns an array of the eight corner points of this 3D rectangle.
// The corners are returned in the order of left to right, back to front,
// bottom to top:
//  left,  back,  bottom
//  right, back,  bottom
//  left,  front, bottom
//  right, front, bottom
//  left,  back,  top
//  right, back,  top
//  left,  front, top
//  right, front, top
func (r Rect3) Corners() [8]Vec3 {
	min := r.Min
	max := r.Max
	return [8]Vec3{
		Vec3{min.X, min.Y, min.Z}, // left, back, bottom
		Vec3{max.X, min.Y, min.Z}, // right, back, bottom
		Vec3{min.X, max.Y, min.Z}, // left, front, bottom
		Vec3{max.X, max.Y, min.Z}, // right, front, bottom

		Vec3{min.X, min.Y, max.Z}, // left, back, top
		Vec3{max.X, min.Y, max.Z}, // right, back, top
		Vec3{min.X, max.Y, max.Z}, // left, front, top
		Vec3{max.X, max.Y, max.Z}, // right, front, top
	}
}

// Area returns the area of this rectangle (the sum of it's sides).
func (r Rect3) Area() float64 {
	s := r.Size()
	return s.X + s.X + s.Y + s.Y + s.Z + s.Z
}

// InSphere reports whether the rectangle r is completely inside the sphere s.
func (r Rect3) InSphere(s Sphere) bool {
	return s.Contains(r.Min) && s.Contains(r.Max)
}

// Rect3Zero is the zero rectangle.
var Rect3Zero = Rect3{}
