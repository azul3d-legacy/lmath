// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
	"math"
)

// Quat represents a four component vector.
type Quat struct {
	W, X, Y, Z float64
}

// String returns an string representation of this vector.
func (a Quat) String() string {
	return fmt.Sprintf("Quat(W=%f, X=%f, Y=%f, Z=%f)", a.W, a.X, a.Y, a.Z)
}

// AlmostEquals tells if a == b using the specified epsilon value.
func (a Quat) AlmostEquals(b Quat, epsilon float64) bool {
	return AlmostEqual(a.W, b.W, epsilon) && AlmostEqual(a.X, b.X, epsilon) && AlmostEqual(a.Y, b.Y, epsilon) && AlmostEqual(a.Z, b.Z, epsilon)
}

// Equals tells if a == b using the default EPSILON value.
func (a Quat) Equals(b Quat) bool {
	return a.AlmostEquals(b, EPSILON)
}

// Add performs a componentwise addition of the two quaternions, returning
// a + b.
func (a Quat) Add(b Quat) Quat {
	return Quat{a.W + a.W, a.X + b.X, a.Y + b.Y, a.Z + b.Z}
}

// AddScalar performs a componentwise scalar addition of a + b.
func (a Quat) AddScalar(b float64) Quat {
	return Quat{a.W + b, a.X + b, a.Y + b, a.Z + b}
}

// Sub performs a componentwise subtraction of the two quaternions, returning
// a - b.
func (a Quat) Sub(b Quat) Quat {
	return Quat{a.W - b.W, a.X - b.X, a.Y - b.Y, a.Z - b.Z}
}

// SubScalar performs a componentwise scalar subtraction of a - b.
func (a Quat) SubScalar(b float64) Quat {
	return Quat{a.W - b, a.X - b, a.Y - b, a.Z - b}
}

// Mul returns the result of the quaternion multiplication a * b.
func (a Quat) Mul(b Quat) Quat {
	return Quat{
		(b.W * a.W) - (b.X * a.X) - (b.Y * a.Y) - (b.Z * a.Z),
		(b.X * a.W) + (b.W * a.X) - (b.Z * a.Y) + (b.Y * a.Z),
		(b.Y * a.W) + (b.Z * a.X) + (b.W * a.Y) - (b.X * a.Z),
		(b.Z * a.W) - (b.Y * a.X) + (b.X * a.Y) + (b.W * a.Z),
	}
}

// MulScalar performs a componentwise scalar multiplication of a * b.
func (a Quat) MulScalar(b float64) Quat {
	return Quat{a.W * b, a.X * b, a.Y * b, a.Z * b}
}

// Div performs a componentwise division of the two quaternions, returning
// a * b.
func (a Quat) Div(b Quat) Quat {
	return Quat{a.W / b.W, a.X / b.X, a.Y / b.Y, a.Z / b.Z}
}

// DivScalar performs a componentwise scalar division of a * b.
func (a Quat) DivScalar(b float64) Quat {
	return Quat{a.W / b, a.X / b, a.Y / b, a.Z / b}
}

// IsNaN tells if any components of this quaternion are not an number.
func (a Quat) IsNaN() bool {
	return math.IsNaN(a.W) || math.IsNaN(a.X) || math.IsNaN(a.Y) || math.IsNaN(a.Z)
}

// Clamp clamps each value in the quaternion to the range of [min, max] and returns
// it.
func (a Quat) Clamp(min, max float64) Quat {
	return Quat{
		Clamp(a.W, min, max),
		Clamp(a.X, min, max),
		Clamp(a.Y, min, max),
		Clamp(a.Z, min, max),
	}
}

// Dot returns the dot product of a and b.
func (a Quat) Dot(b Quat) float64 {
	return a.W*b.W + a.X*b.X + a.Y*b.Y + a.Z*b.Z
}

// Inverse returns the inverse of this quaternion.
func (a Quat) Inverse() Quat {
	return Quat{
		-a.W,
		a.X,
		a.Y,
		a.Z,
	}
}

// LengthSq returns the magnitude squared of this quaternion, useful for comparing
// distances.
func (a Quat) LengthSq() float64 {
	return a.W*a.W + a.X*a.X + a.Y*a.Y + a.Z*a.Z
}

// Length returns the magnitude of this quaternion. To avoid a sqrt call when
// strictly comparing distances, LengthSq can be used instead.
func (a Quat) Length() float64 {
	return math.Sqrt(a.W*a.W + a.X*a.X + a.Y*a.Y + a.Z*a.Z)
}

// Normalized returns the normalized (i.e. length/magnitude == 1) quaternion of
// a.
// The quaternion a must be non-zero or else division by zero may occur.
func (a Quat) Normalized() Quat {
	length := math.Sqrt(a.W*a.W + a.X*a.X + a.Y*a.Y + a.Z*a.Z)
	return Quat{
		a.W / length,
		a.X / length,
		a.Y / length,
		a.Z / length,
	}
}

// Min returns a quaternion representing the smallest components of both the
// quaternions.
func (a Quat) Min(b Quat) Quat {
	var r Quat
	if a.W < b.W {
		r.W = a.W
	} else {
		r.W = b.W
	}
	if a.X < b.X {
		r.X = a.X
	} else {
		r.X = b.X
	}
	if a.Y < b.Y {
		r.Y = a.Y
	} else {
		r.Y = b.Y
	}
	if a.Z < b.Z {
		r.Z = a.Z
	} else {
		r.Z = b.Z
	}
	return r
}

// Max returns a quaternion representing the largest components of both the
// quaternions.
func (a Quat) Max(b Quat) Quat {
	var r Quat
	if a.W > b.W {
		r.W = a.W
	} else {
		r.W = b.W
	}
	if a.X > b.X {
		r.X = a.X
	} else {
		r.X = b.X
	}
	if a.Y > b.Y {
		r.Y = a.Y
	} else {
		r.Y = b.Y
	}
	if a.Z > b.Z {
		r.Z = a.Z
	} else {
		r.Z = b.Z
	}
	return r
}

// Lerp returns a quaternion representing the linear interpolation between the
// quaternions a and b. The t parameter is the amount to interpolate
// (0.0 - 1.0) between the quaternions.
func (a Quat) Lerp(b Quat, t float64) Quat {
	return a.Mul(b.MulScalar(t))
}

// Conjugate calculates and returns the conjugate of this quaternion.
func (a Quat) Conjugate() Quat {
	return Quat{
		a.W,
		-a.X,
		-a.Y,
		-a.Z,
	}
}

// Hpr extracts the equivilent heading, pitch, and roll euler angles in radians
// from the quaternion for the given coordinate system.
func (a Quat) Hpr(cs CoordSys) (hpr Vec3) {
	mat := a.ExtractToMat3()
	_, _, hpr = mat.Decompose(cs)
	return
}

// Up returns the orientation represented by this quaternion expressed as an up
// vector.
func (a Quat) Up(cs CoordSys) Vec3 {
	return a.TransformVec3(cs.Up())
}

// Right returns the orientation represented by this quaternion expressed as an
// right vector.
func (a Quat) Right(cs CoordSys) Vec3 {
	return a.TransformVec3(cs.Right())
}

// Forward returns the orientation represented by this quaternion expressed as
// an forward vector.
func (a Quat) Forward(cs CoordSys) Vec3 {
	return a.TransformVec3(cs.Forward())
}

// Angle returns the angle between the orientation represented by this
// quaternion and the orientation represent by the quaternion b in radians.
func (a Quat) AngleQuat(cs CoordSys, b Quat) float64 {
	return a.Forward(cs).Angle(b.Forward(cs))
}

// Axis returns the axis of the rotation represented by the quaternion. The
// returned vector is not normalized.
func (a Quat) Axis() Vec3 {
	return Vec3{a.X, a.Y, a.Z}
}

// Angle returns the rotation represented by the quaternion as an angle about
// an arbitrary axis (returned by the Axis() function), the returned value is
// expressed in radian units counterclockwise about the axis.
func (a Quat) Angle() float64 {
	a = a.Normalized()
	return math.Acos(a.W) * 2
}

// TransformVec3 transforms the 3-component vector by the specified quaternion
// rotation and returns the result.
func (a Quat) TransformVec3(v Vec3) Vec3 {
	vecQuat := Quat{0.0, v.X, v.Y, v.Z}
	conjugate := a.Conjugate()
	vecQuat = conjugate.Mul(vecQuat).Mul(a)
	return Vec3{vecQuat.X, vecQuat.Y, vecQuat.Z}
}

// TransformVec4 transforms the 4-component vector by the specified quaternion
// rotation and returns the result.
func (a Quat) TransformVec4(v Vec4) Vec4 {
	vecQuat := Quat{v.W, v.X, v.Y, v.Z}
	conjugate := a.Conjugate()
	vecQuat = conjugate.Mul(vecQuat).Mul(a)
	return vecQuat.Vec4()
}

// Vec4 converts the quaternion into a four-component vector; Short-hand for:
//  Vec4{a.W, a.X, a.Y, a.Z}
func (a Quat) Vec4() Vec4 {
	return Vec4{a.W, a.X, a.Y, a.Z}
}

// ExtractToMat3 extracts the quaternion into a three-component matrix and
// returns it.
func (a Quat) ExtractToMat3() Mat3 {
	n := a.Dot(a)

	var s float64
	if !Equal(n, 0) {
		s = 2.0 / n
	}

	xs := a.X * s
	ys := a.Y * s
	zs := a.Z * s
	wx := a.W * xs
	wy := a.W * ys
	wz := a.W * zs
	xx := a.X * xs
	xy := a.X * ys
	xz := a.X * zs
	yy := a.Y * ys
	yz := a.Y * zs
	zz := a.Z * zs

	return Matrix3(
		1.0-(yy+zz), xy+wz, xz-wy,
		xy-wz, 1.0-(xx+zz), yz+wx,
		xz+wy, yz-wx, 1.0-(xx+yy),
	)
}

// ExtractToMat4 extracts the quaternion into a four-component matrix and
// returns it.
func (a Quat) ExtractToMat4() Mat4 {
	n := a.Dot(a)

	var s float64
	if !Equal(n, 0) {
		s = 2.0 / n
	}

	xs := a.X * s
	ys := a.Y * s
	zs := a.Z * s
	wx := a.W * xs
	wy := a.W * ys
	wz := a.W * zs
	xx := a.X * xs
	xy := a.X * ys
	xz := a.X * zs
	yy := a.Y * ys
	yz := a.Y * zs
	zz := a.Z * zs

	return Matrix4(
		1.0-(yy+zz), xy+wz, xz-wy, 0,
		xy-wz, 1.0-(xx+zz), yz+wx, 0,
		xz+wy, yz-wx, 1.0-(xx+yy), 0,
		0, 0, 0, 1,
	)
}

// QuatFromMat3 returns a quaternion rotation according to the rotation
// represented by the matrix.
func QuatFromMat3(m Mat3) Quat {
	var a Quat
	m00 := m[0][0]
	m10 := m[1][0]
	m20 := m[2][0]
	m01 := m[0][1]
	m11 := m[1][1]
	m21 := m[2][1]
	m02 := m[0][2]
	m12 := m[1][2]
	m22 := m[2][2]

	trace := m00 + m11 + m22

	if trace > 0 {
		// The easy case.
		s := math.Sqrt(trace + 1.0)
		a.W = s * 0.5
		s = 0.5 / s
		a.X = (m12 - m21) * s
		a.Y = (m20 - m02) * s
		a.Z = (m01 - m10) * s
	} else {
		// The harder case.  First, figure out which column to take as
		// root.  This will be the column with the largest value.
		//
		// It is tempting to try to compare the absolute values of the
		// diagonal values in the code below, instead of their normal,
		// signed values.  Don't do it.  We are actually maximizing the
		// value of S, which must always be positive, and is therefore
		// based on the diagonal whose actual value--not absolute
		// value--is greater than those of the other two.
		//
		// We already know that m00 + m11 + m22 <= 0 (because we are here
		// in the harder case).
		if m00 > m11 && m00 > m22 {
			// m00 is larger than m11 and m22.
			s := 1.0 + m00 - (m11 + m22)
			s = math.Sqrt(s)
			a.X = s * 0.5
			s = 0.5 / s
			a.Y = (m01 + m10) * s
			a.Z = (m02 + m20) * s
			a.W = (m12 - m21) * s

		} else if m11 > m22 {
			// m11 is larger than m00 and m22.
			s := 1.0 + m11 - (m22 + m00)
			s = math.Sqrt(s)
			a.Y = s * 0.5
			s = 0.5 / s
			a.Z = (m12 + m21) * s
			a.X = (m10 + m01) * s
			a.W = (m20 - m02) * s

		} else {
			// m22 is larger than m00 and m11.
			s := 1.0 + m22 - (m00 + m11)
			s = math.Sqrt(s)
			a.Z = s * 0.5
			s = 0.5 / s
			a.X = (m20 + m02) * s
			a.Y = (m21 + m12) * s
			a.W = (m01 - m10) * s
		}
	}
	return a
}

// QuatFromHpr returns a quaternion equivilent to the heading, pithc, and roll
// Euler angles in radians for the given coordinate system.
func QuatFromHpr(hpr Vec3, cs CoordSys) Quat {
	v := cs.Up()
	n := hpr.X * 0.5
	s := math.Sin(n)
	c := math.Cos(n)
	quatHeading := Quat{c, v.X * s, v.Y * s, v.Z * s}

	v = cs.Right()
	n = hpr.Y * 0.5
	s = math.Sin(n)
	c = math.Cos(n)
	quatPitch := Quat{c, v.X * s, v.Y * s, v.Z * s}

	v = cs.Forward()
	n = hpr.Z * 0.5
	s = math.Sin(n)
	c = math.Cos(n)
	quatRoll := Quat{c, v.X * s, v.Y * s, v.Z * s}

	if cs.RightHanded() {
		return quatRoll.Mul(quatPitch).Mul(quatHeading)
	}
	return quatHeading.Mul(quatPitch).Mul(quatRoll).Inverse()
}

// QuatFromAxisAngle returns a quaternion rotation which represents the given
// angle of rotation in radians about the given axis.
func QuatFromAxisAngle(axis Vec3, angle float64) Quat {
	axis, _ = axis.Normalized()
	sinHalfAngle := math.Sin(angle * 0.5)
	return Quat{
		math.Cos(angle * 0.5),
		axis.X * sinHalfAngle,
		axis.Y * sinHalfAngle,
		axis.Z * sinHalfAngle,
	}
}

var (
	QuatIdentity = Quat{1, 0, 0, 0}
	QuatZero     = Quat{0, 0, 0, 0}
)
