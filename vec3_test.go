// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestVec3Add(t *testing.T) {
	a := Vec3{1, 3, 3}
	b := Vec3{1, 3, 3}
	if !a.Add(b).Equals(Vec3{2, 6, 6}) {
		t.Fail()
	}

	if !a.AddScalar(2).Equals(Vec3{3, 5, 5}) {
		t.Fail()
	}
}

func TestVec3Sub(t *testing.T) {
	a := Vec3{1, 3, 3}
	b := Vec3{1, 3, 3}
	if !a.Sub(b).Equals(Vec3{0, 0, 0}) {
		t.Fail()
	}

	if !a.SubScalar(2).Equals(Vec3{-1, 1, 1}) {
		t.Fail()
	}
}

func TestVec3Mul(t *testing.T) {
	a := Vec3{1, 3, 3}
	b := Vec3{1, 3, 3}
	if !a.Mul(b).Equals(Vec3{1, 9, 9}) {
		t.Fail()
	}

	if !a.MulScalar(2).Equals(Vec3{2, 6, 6}) {
		t.Fail()
	}
}

func TestVec3Div(t *testing.T) {
	a := Vec3{1, 3, 3}
	b := Vec3{1, 3, 3}
	if !a.Div(b).Equals(Vec3{1, 1, 1}) {
		t.Fail()
	}

	if !a.DivScalar(2).Equals(Vec3{0.5, 3.0 / 2.0, 3.0 / 2.0}) {
		t.Fail()
	}
}

func BenchmarkVec3Equals(b *testing.B) {
	x := Vec3{1, 3, 3}
	y := Vec3{1.33, 3.33, 3.33}
	for n := 0; n < b.N; n++ {
		x.Equals(y)
	}
}
