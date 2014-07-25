// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestMat4Row(t *testing.T) {
	x := Matrix4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	)
	if !x.Row(0).Equals(Vec4{1, 0, 0, 0}) {
		t.Fail()
	}
	if !x.Row(1).Equals(Vec4{0, 1, 0, 0}) {
		t.Fail()
	}
	if !x.Row(2).Equals(Vec4{0, 0, 1, 0}) {
		t.Fail()
	}
	if !x.Row(3).Equals(Vec4{0, 0, 0, 1}) {
		t.Fail()
	}

	x = x.SetRow(0, Vec4{0, 1, 0, 0})
	x = x.SetRow(1, Vec4{0, 0, 1, 0})
	x = x.SetRow(2, Vec4{1, 0, 0, 0})
	x = x.SetRow(3, Vec4{0, 0, 0, 1})
	supposed := Matrix4(
		0, 1, 0, 0,
		0, 0, 1, 0,
		1, 0, 0, 0,
		0, 0, 0, 1,
	)
	if !x.Equals(supposed) {
		t.Fail()
	}
}
