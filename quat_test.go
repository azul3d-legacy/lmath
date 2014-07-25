// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestQuatToMat3(t *testing.T) {
	before := Mat3FromAxisAngle(Vec3{0, 1, 0}, Radians(45), CoordSysZUpRight)
	quat := QuatFromMat3(before)
	after := quat.ExtractToMat3()

	if !before.Equals(after) {
		t.Log(before)
		t.Log(after)
		t.Fail()
	}
}

func TestQuatToEuler(t *testing.T) {
	hpr := Vec3{0, 45, 0}.XyzToHpr().Radians()
	before := QuatFromHpr(hpr, CoordSysZUpRight)
	after := before.Hpr(CoordSysZUpRight)
	if !hpr.Equals(after) {
		t.Log(hpr)
		t.Log(after)
		t.Fail()
	}
}
