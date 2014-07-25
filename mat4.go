// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"fmt"
	"math"
)

// Mat4 represents an 4x4 matrix, indices are in m[row][column] order.
type Mat4 [4][4]float64

// String returns an string representation of this matrix.
func (a Mat4) String() string {
	return fmt.Sprintf("Mat4(\n    %f, %f, %f, %f\n    %f, %f, %f, %f\n    %f, %f, %f, %f\n    %f, %f, %f, %f\n)", a[0][0], a[0][1], a[0][2], a[0][3], a[1][0], a[1][1], a[1][2], a[1][3], a[2][0], a[2][1], a[2][2], a[2][3], a[3][0], a[3][1], a[3][2], a[3][3])
}

// SetRow sets the values in the specified matrix row to the values in the
// specified four-component vector.
// The row parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) SetRow(row int, values Vec4) Mat4 {
	a[row][0] = values.X
	a[row][1] = values.Y
	a[row][2] = values.Z
	a[row][3] = values.W
	return a
}

// Row returns the values in the specified matrix row as an three-component
// vector.
// The row parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) Row(row int) Vec4 {
	return Vec4{a[row][0], a[row][1], a[row][2], a[row][3]}
}

// SetRowVec3 sets the values in the specified matrix row to the values in the
// specified three-component vector, leaving the fourth element in the row
// untouched.
// The row parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) SetRowVec3(row int, values Vec3) Mat4 {
	a[row][0] = values.X
	a[row][1] = values.Y
	a[row][2] = values.Z
	return a
}

// RowVec3 returns the values in the specified matrix row as an three-component
// vector, the fourth element of the row is ignored.
// The row parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) RowVec3(row int) Vec3 {
	return Vec3{a[row][0], a[row][1], a[row][2]}
}

// SetCol sets the values in the specified matrix column to the values in the
// specified four-component vector.
// The column parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) SetCol(column int, values Vec4) Mat4 {
	a[0][column] = values.X
	a[1][column] = values.Y
	a[2][column] = values.Z
	a[3][column] = values.W
	return a
}

// Col returns the values in the specified matrix column as an four-component
// vector.
// The column parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) Col(column int) Vec4 {
	return Vec4{a[0][column], a[1][column], a[2][column], a[3][column]}
}

// SetColVec3 sets the values in the specified matrix column to the values in
// the specified three-component vector, leaving the fourth element in the
// column untouched.
// The column parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) SetColVec3(column int, values Vec3) Mat4 {
	a[0][column] = values.X
	a[1][column] = values.Y
	a[2][column] = values.Z
	return a
}

// ColVec3 returns the values in the specified matrix column as an
// three-component vector, ignoring the fourth element of the column.
// The column parameter must be 0, 1, 2, or 3 or else an panic will occur.
func (a Mat4) ColVec3(column int) Vec3 {
	return Vec3{a[0][column], a[1][column], a[2][column]}
}

// AlmostEquals tells whether a is memberwise equal to b using the specified
// epsilon value.
func (a Mat4) AlmostEquals(b Mat4, epsilon float64) bool {
	return AlmostEqual(a[0][0], b[0][0], epsilon) &&
		AlmostEqual(a[0][1], b[0][1], epsilon) &&
		AlmostEqual(a[0][2], b[0][2], epsilon) &&
		AlmostEqual(a[0][3], b[0][3], epsilon) &&

		AlmostEqual(a[1][0], b[1][0], epsilon) &&
		AlmostEqual(a[1][1], b[1][1], epsilon) &&
		AlmostEqual(a[1][2], b[1][2], epsilon) &&
		AlmostEqual(a[1][3], b[1][3], epsilon) &&

		AlmostEqual(a[2][0], b[2][0], epsilon) &&
		AlmostEqual(a[2][1], b[2][1], epsilon) &&
		AlmostEqual(a[2][2], b[2][2], epsilon) &&
		AlmostEqual(a[2][3], b[2][3], epsilon) &&

		AlmostEqual(a[3][0], b[3][0], epsilon) &&
		AlmostEqual(a[3][1], b[3][1], epsilon) &&
		AlmostEqual(a[3][2], b[3][2], epsilon) &&
		AlmostEqual(a[3][3], b[3][3], epsilon)
}

// Equals tells whether a is memberwise equal to b using the default EPSILON
// value.
func (a Mat4) Equals(b Mat4) bool {
	return a.AlmostEquals(b, EPSILON)
}

// AddScalar performs in-place memberwise addition a + s and returns the
// result.
func (a Mat4) AddScalar(s float64) Mat4 {
	return Matrix4(
		a[0][0]+s, a[0][1]+s, a[0][2]+s, a[0][3]+s,
		a[1][0]+s, a[1][1]+s, a[1][2]+s, a[1][3]+s,
		a[2][0]+s, a[2][1]+s, a[2][2]+s, a[2][3]+s,
		a[3][0]+s, a[3][1]+s, a[3][2]+s, a[3][3]+s,
	)
}

// SubScalar performs in-place memberwise subtraction a - s and returns the
// result.
func (a Mat4) SubScalar(s float64) Mat4 {
	return Matrix4(
		a[0][0]-s, a[0][1]-s, a[0][2]-s, a[0][3]-s,
		a[1][0]-s, a[1][1]-s, a[1][2]-s, a[1][3]-s,
		a[2][0]-s, a[2][1]-s, a[2][2]-s, a[2][3]-s,
		a[3][0]-s, a[3][1]-s, a[3][2]-s, a[3][3]-s,
	)
}

// MulScalar performs in-place memberwise multiplication a * s and returns the
// result.
func (a Mat4) MulScalar(s float64) Mat4 {
	return Matrix4(
		a[0][0]*s, a[0][1]*s, a[0][2]*s, a[0][3]*s,
		a[1][0]*s, a[1][1]*s, a[1][2]*s, a[1][3]*s,
		a[2][0]*s, a[2][1]*s, a[2][2]*s, a[2][3]*s,
		a[3][0]*s, a[3][1]*s, a[3][2]*s, a[3][3]*s,
	)
}

// DivScalar performs in-place memberwise division a / s and returns the
// result.
func (a Mat4) DivScalar(s float64) Mat4 {
	return Matrix4(
		a[0][0]/s, a[0][1]/s, a[0][2]/s, a[0][3]/s,
		a[1][0]/s, a[1][1]/s, a[1][2]/s, a[1][3]/s,
		a[2][0]/s, a[2][1]/s, a[2][2]/s, a[2][3]/s,
		a[3][0]/s, a[3][1]/s, a[3][2]/s, a[3][3]/s,
	)
}

// Add performs in-place memberwise addition a + b and returns the result.
func (a Mat4) Add(b Mat4) Mat4 {
	return Matrix4(
		a[0][0]+b[0][0], a[0][1]+b[0][1], a[0][2]+b[0][2], a[0][3]+b[0][3],
		a[1][0]+b[1][0], a[1][1]+b[1][1], a[1][2]+b[1][2], a[1][3]+b[1][3],
		a[2][0]+b[2][0], a[2][1]+b[2][1], a[2][2]+b[2][2], a[2][3]+b[2][3],
		a[3][0]+b[3][0], a[3][1]+b[3][1], a[3][2]+b[3][2], a[3][3]+b[3][3],
	)
}

// Sub performs in-place memberwise subtraction a - b and returns the result.
func (a Mat4) Sub(b Mat4) Mat4 {
	return Matrix4(
		a[0][0]-b[0][0], a[0][1]-b[0][1], a[0][2]-b[0][2], a[0][3]-b[0][3],
		a[1][0]-b[1][0], a[1][1]-b[1][1], a[1][2]-b[1][2], a[1][3]-b[1][3],
		a[2][0]-b[2][0], a[2][1]-b[2][1], a[2][2]-b[2][2], a[2][3]-b[2][3],
		a[3][0]-b[3][0], a[3][1]-b[3][1], a[3][2]-b[3][2], a[3][3]-b[3][3],
	)
}

// Mul performs matrix multiplication a * b and returns the result.
func (a Mat4) Mul(b Mat4) Mat4 {
	var out Mat4
	out[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0] + a[0][3]*b[3][0]
	out[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1] + a[0][3]*b[3][1]
	out[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2] + a[0][3]*b[3][2]
	out[0][3] = a[0][0]*b[0][3] + a[0][1]*b[1][3] + a[0][2]*b[2][3] + a[0][3]*b[3][3]

	out[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0] + a[1][3]*b[3][0]
	out[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1] + a[1][3]*b[3][1]
	out[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2] + a[1][3]*b[3][2]
	out[1][3] = a[1][0]*b[0][3] + a[1][1]*b[1][3] + a[1][2]*b[2][3] + a[1][3]*b[3][3]

	out[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0] + a[2][3]*b[3][0]
	out[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1] + a[2][3]*b[3][1]
	out[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2] + a[2][3]*b[3][2]
	out[2][3] = a[2][0]*b[0][3] + a[2][1]*b[1][3] + a[2][2]*b[2][3] + a[2][3]*b[3][3]

	out[3][0] = a[3][0]*b[0][0] + a[3][1]*b[1][0] + a[3][2]*b[2][0] + a[3][3]*b[3][0]
	out[3][1] = a[3][0]*b[0][1] + a[3][1]*b[1][1] + a[3][2]*b[2][1] + a[3][3]*b[3][1]
	out[3][2] = a[3][0]*b[0][2] + a[3][1]*b[1][2] + a[3][2]*b[2][2] + a[3][3]*b[3][2]
	out[3][3] = a[3][0]*b[0][3] + a[3][1]*b[1][3] + a[3][2]*b[2][3] + a[3][3]*b[3][3]
	return out
}

// Transposed returns the transposed version of matrix a.
func (a Mat4) Transposed() Mat4 {
	a01 := a[0][1]
	a[0][1] = a[1][0]
	a[1][0] = a01

	a02 := a[0][2]
	a[0][2] = a[2][0]
	a[2][0] = a02

	a03 := a[0][3]
	a[0][3] = a[3][0]
	a[3][0] = a03

	a12 := a[1][2]
	a[1][2] = a[2][1]
	a[2][1] = a12

	a13 := a[1][3]
	a[1][3] = a[3][1]
	a[3][1] = a13

	a23 := a[2][3]
	a[2][3] = a[3][2]
	a[3][2] = a23
	return a
}

// SetUpperMat3 sets the upper-left 3x3 matrix to the specified one and returns
// the new 4x4 matrix.
func (a Mat4) SetUpperMat3(b Mat3) Mat4 {
	a[0][0] = b[0][0]
	a[0][1] = b[0][1]
	a[0][2] = b[0][2]

	a[1][0] = b[1][0]
	a[1][1] = b[1][1]
	a[1][2] = b[1][2]

	a[2][0] = b[2][0]
	a[2][1] = b[2][1]
	a[2][2] = b[2][2]
	return a
}

// UpperMat3 returns the upper-left 3x3 matrix.
func (a Mat4) UpperMat3() Mat3 {
	return Matrix3(
		a[0][0], a[0][1], a[0][2],
		a[1][0], a[1][1], a[1][2],
		a[2][0], a[2][1], a[2][2],
	)
}

// MulQuat multiplies the matrix by the quaternion a * b and returns the
// result.
func (a Mat4) MulQuat(b Quat) Mat4 {
	quatMat := b.ExtractToMat4()

	// Preserve the homogeneous coords and the translate
	row3 := a.Row(3)
	col3 := a.Col(3)

	quatMat = a.Mul(quatMat)
	quatMat = quatMat.SetRow(3, row3)
	quatMat = quatMat.SetCol(3, col3)
	return quatMat
}

// AffineInverse returns the inverse of the affine transformation matrix and
// returns the result.
// Returned is Mat4Identity, ok=false if the matrix was singular and could not
// be inverted,
func (a Mat4) AffineInverse() (out Mat4, ok bool) {
	rot, ok := a.UpperMat3().Inverse()
	if !ok {
		return Mat4Identity, false
	}

	out = Mat4Identity.SetUpperMat3(rot)
	out[3][0] = -(a[3][0]*out[0][0] + a[3][1]*out[1][0] + a[3][2]*out[2][0])
	out[3][1] = -(a[3][0]*out[0][1] + a[3][1]*out[1][1] + a[3][2]*out[2][1])
	out[3][2] = -(a[3][0]*out[0][2] + a[3][1]*out[1][2] + a[3][2]*out[2][2])
	return out, true
}

func (m *Mat4) decomposeMat(index *[4]int) bool {
	var vv [4]float64

	for i := 0; i < 4; i++ {
		big := 0.0
		for j := 0; j < 4; j++ {
			temp := math.Abs(m[i][j])
			if temp > big {
				big = temp
			}
		}

		// We throw the value out only if it's smaller than our "small"
		// threshold squared.  This helps reduce overly-sensitive
		// rejections.
		if Equal(big, EPSILON*EPSILON) {
			return false
		}
		vv[i] = 1.0 / big
	}

	for j := 0; j < 4; j++ {
		for i := 0; i < j; i++ {
			sum := m[i][j]
			for k := 0; k < i; k++ {
				sum -= m[i][k] * m[k][j]
			}
			m[i][j] = sum
		}

		big := 0.0
		imax := -1

		for i := j; i < 4; i++ {
			sum := m[i][j]
			for k := 0; k < j; k++ {
				sum -= m[i][k] * m[k][j]
			}

			m[i][j] = sum

			dum := vv[i] * math.Abs(sum)

			if dum >= big {
				big = dum
				imax = i
			}
		}
		if j != imax {
			for k := 0; k < 4; k++ {
				dum := m[imax][k]
				m[imax][k] = m[j][k]
				m[j][k] = dum
			}
			vv[imax] = vv[j]
		}
		index[j] = imax

		if Equal(m[j][j], 0) {
			m[j][j] = EPSILON
		}

		if j != 4-1 {
			dum := 1.0 / m[j][j]

			for i := j + 1; i < 4; i++ {
				m[i][j] *= dum
			}
		}
	}
	return true
}

func (m *Mat4) backSubMat(index *[4]int, inv *Mat4, row int) bool {
	var i int

	ii := -1

	for i = 0; i < 4; i++ {
		ip := index[i]
		sum := inv[row][ip]
		inv[row][ip] = inv[row][i]
		if ii >= 0 {
			for j := ii; j <= i-1; j++ {
				sum -= m[i][j] * inv[row][j]
			}
		} else if sum != 0 {
			ii = i
		}

		inv[row][i] = sum
	}

	for i = 4 - 1; i >= 0; i-- {
		sum := inv[row][i]
		for j := i + 1; j < 4; j++ {
			sum -= m[i][j] * inv[row][j]
		}
		inv[row][i] = sum / m[i][i]
	}

	return true
}

// Inverse returns the inverse of the matrix, as a fully general operation and
// makes no requirements about the type of transform represented by the matrix.
// Returned is Mat4Identity, ok=false if the matrix was singular and could not
// be inverted.
func (a Mat4) Inverse() (result Mat4, ok bool) {
	if Equal(a[0][3], 0) && Equal(a[1][3], 0) && Equal(a[2][3], 0) && Equal(a[3][3], 0) {
		return a.AffineInverse()
	}

	result = a
	out := &result
	var index [4]int

	// Other order?
	if !out.decomposeMat(&index) {
		return Mat4Identity, false
	}

	inv := Mat4Identity
	for row := 0; row < 4; row++ {
		out.backSubMat(&index, &inv, row)
	}

	return inv.Transposed(), true
}

// SetTranslation sets the translation components of the matrix to the given
// vector and returns the result.
func (a Mat4) SetTranslation(t Vec3) Mat4 {
	a[3][0] = t.X
	a[3][1] = t.Y
	a[3][2] = t.Z
	return a
}

// Translation returns the translation components of the matrix as a vector.
func (a Mat4) Translation() Vec3 {
	return Vec3{a[3][0], a[3][1], a[3][2]}
}

// IsNan tells if any components of this matrix are not an number.
func (a Mat4) IsNaN() bool {
	return math.IsNaN(a[0][0]) || math.IsNaN(a[0][1]) || math.IsNaN(a[0][2]) || math.IsNaN(a[0][3]) ||
		math.IsNaN(a[1][0]) || math.IsNaN(a[1][1]) || math.IsNaN(a[1][2]) || math.IsNaN(a[1][3]) ||
		math.IsNaN(a[2][0]) || math.IsNaN(a[2][1]) || math.IsNaN(a[2][2]) || math.IsNaN(a[2][3]) ||
		math.IsNaN(a[3][0]) || math.IsNaN(a[3][1]) || math.IsNaN(a[3][2]) || math.IsNaN(a[3][3])
}

// Project returns a 2D point in the range -1 to +1 given a 3D point also in
// the frustum matrix a's coordinate space.
//
// If ok=false is returned then the point is outside of the frustum matrix a,
// and the returned point may not be meaningful.
func (a Mat4) Project(p3 Vec3) (p2 Vec2, ok bool) {
	p4 := Vec4{p3.X, p3.Y, p3.Z, 1.0}
	p4 = p4.Transform(a)
	if p4.W == 0 {
		p2 = Vec2Zero
		ok = false
		return
	}

	recipW := 1.0 / p4.W
	p2 = Vec2{p4.X * recipW, p4.Y * recipW}

	xValid := (p2.X >= -1) && (p2.X <= 1)
	yValid := (p2.Y >= -1) && (p2.Y <= 1)
	ok = (p4.W > 0) && xValid && yValid
	return
}

// Mat4FromAxisAngle returns a rotation matrix that will rotate by the given
// angle in radians counterclockwise about the indicated axis.
func Mat4FromAxisAngle(axis Vec3, angle float64, cs CoordSys) Mat4 {
	var a Mat4
	axis, _ = axis.Normalized()

	// In a left-handed coordinate system, counterclockwise is the
	// other direction.
	if cs.LeftHanded() {
		angle = -angle
	}

	s := math.Sin(angle)
	c := math.Cos(angle)
	t := 1.0 - c

	t0 := t * axis.X
	t1 := t * axis.Y
	t2 := t * axis.Z
	s0 := s * axis.X
	s1 := s * axis.Y
	s2 := s * axis.Z

	a[0][0] = t0*axis.X + c
	a[0][1] = t0*axis.Y + s2
	a[0][2] = t0*axis.Z - s1

	a[1][0] = t1*axis.X - s2
	a[1][1] = t1*axis.Y + c
	a[1][2] = t1*axis.Z + s0

	a[2][0] = t2*axis.X + s1
	a[2][1] = t2*axis.Y - s0
	a[2][2] = t2*axis.Z + c

	a[0][3] = 0
	a[1][3] = 0
	a[2][3] = 0

	a[3][0] = 0
	a[3][1] = 0
	a[3][2] = 0
	a[3][3] = 1
	return a
}

// Mat4FromScaleShear returns a matrix that will apply the given scaling and
// shearing values along their respective axis in the specified coordinate
// system.
// A panic will occur if the coordinate system is invalid.
func Mat4FromScaleShear(scale, shear Vec3, cs CoordSys) Mat4 {
	m3 := Mat3FromScaleShear(scale, shear, cs)
	return Mat4Identity.SetUpperMat3(m3)
}

// Mat4FromTranslation returns a matrix that will apply the given translation
// vector.
func Mat4FromTranslation(translation Vec3) Mat4 {
	return Matrix4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		translation.X, translation.Y, translation.Z, 1,
	)
}

// Mat4FromScale returns a matrix that will apply the given scaling values
// along their respective axis.
func Mat4FromScale(scale Vec3) Mat4 {
	return Matrix4(
		scale.X, 0, 0, 0,
		0, scale.Y, 0, 0,
		0, 0, scale.Z, 0,
		0, 0, 0, 1,
	)
}

// Mat4FromFrustum returns a matrix that represents the given frustum given the
// specified bounds.
func Mat4FromFrustum(left, right, bottom, top, near, far float64) Mat4 {
	normalWidth := 1.0 / (right - left)
	normalHeight := 1.0 / (top - bottom)
	normalDepth := 1.0 / (near - far)

	return Matrix4(
		(near*2.0)*normalWidth,
		0,
		0,
		0,

		0,
		(near*2.0)*normalHeight,
		0,
		0,

		(right+left)*normalWidth,
		(top+bottom)*normalHeight,
		(far+near)*normalDepth,
		-1,

		0,
		0,
		(far*near*2.0)*normalDepth,
		0,
	)
}

// Mat4Perspective returns a matrix that represents a perspective viewing
// frustum given the specified field of view, aspect ratio, and near/far
// values.
func Mat4Perspective(fovY, aspectRatio, near, far float64) Mat4 {
	fH := math.Tan(fovY/360*math.Pi) * near
	fW := fH * aspectRatio
	return Mat4FromFrustum(-fW, fW, -fH, fH, near, far)
}

// Mat4Ortho returns a orthographic projection matrix given the specified
// frustum bounds.
// See: http://en.wikipedia.org/wiki/Orthographic_projection_(geometry)
func Mat4Ortho(left, right, bottom, top, near, far float64) Mat4 {
	return Matrix4(
		2.0/(right-left),
		0,
		0,
		0,

		0,
		2.0/(top-bottom),
		0,
		0,

		0,
		0,
		-2.0/(far-near),
		0,

		-(right+left)/(right-left),
		-(top+bottom)/(top-bottom),
		-(far+near)/(far-near),
		1,
	)
}

// Mat4UnOrtho returns a orthographic unprojection matrix given the specified
// frustum bounds.
// See: http://en.wikipedia.org/wiki/Orthographic_projection_(geometry)
func Mat4UnOrtho(left, right, bottom, top, near, far float64) Mat4 {
	return Matrix4(
		(right-left)/2.0,
		0,
		0,
		0,

		0,
		(top-bottom)/2.0,
		0,
		0,

		0,
		0,
		(far-near)/-2.0,
		0,

		(left+right)/2,
		(top+bottom)/2,
		(far+near)/-2,
		1,
	)
}

// Matrix4 returns an new *Mat4 given the specified matrix components.
func Matrix4(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33 float64) Mat4 {
	return Mat4{
		{m00, m01, m02, m03},
		{m10, m11, m12, m13},
		{m20, m21, m22, m23},
		{m30, m31, m32, m33},
	}
}

var (
	Mat4Identity = Matrix4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	)

	Mat4Zeros = Matrix4(
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0,
	)

	Mat4Ones = Matrix4(
		1, 1, 1, 1,
		1, 1, 1, 1,
		1, 1, 1, 1,
		1, 1, 1, 1,
	)

	Mat4YToZUp = Matrix4(
		1, 0, 0, 0,
		0, 0, 1, 0,
		0, -1, 0, 0,
		0, 0, 0, 1,
	)

	Mat4ZToYUp = Matrix4(
		1, 0, 0, 0,
		0, 0, -1, 0,
		0, 1, 0, 0,
		0, 0, 0, 1,
	)

	Mat4FlipY = Matrix4(
		1, 0, 0, 0,
		0, -1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	)

	Mat4FlipZ = Matrix4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, -1, 0,
		0, 0, 0, 1,
	)

	Mat4LZToRY = Mat4FlipY.Mul(Mat4ZToYUp)
	Mat4LYToRZ = Mat4FlipZ.Mul(Mat4YToZUp)
)
