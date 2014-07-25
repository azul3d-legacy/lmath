// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestMat3Row(t *testing.T) {
	x := Matrix3(
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	)
	if !x.Row(0).Equals(Vec3{1, 0, 0}) {
		t.Fail()
	}
	if !x.Row(1).Equals(Vec3{0, 1, 0}) {
		t.Fail()
	}
	if !x.Row(2).Equals(Vec3{0, 0, 1}) {
		t.Fail()
	}

	x = x.SetRow(0, Vec3{0, 1, 0})
	x = x.SetRow(1, Vec3{0, 0, 1})
	x = x.SetRow(2, Vec3{1, 0, 0})
	supposed := Matrix3(
		0, 1, 0,
		0, 0, 1,
		1, 0, 0,
	)
	if !x.Equals(supposed) {
		t.Fail()
	}
}

func TestMat3ComposeDecompose(t *testing.T) {
	cs := CoordSysZUpRight
	scale := Vec3{1, 2, 3}
	shear := Vec3{4, 5, 6}
	hpr := Vec3{5, 4, 3}
	m := Mat3Compose(scale, shear, hpr, cs)
	dScale, dShear, dHpr := m.Decompose(cs)
	t.Log(scale, dScale)
	t.Log(shear, dShear)
	t.Log(hpr, dHpr)
	if !scale.AlmostEquals(dScale, 0.1) {
		t.Fail()
	}
	if !shear.AlmostEquals(dShear, 0.1) {
		t.Fail()
	}
}
