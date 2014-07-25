// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestVec4Add(t *testing.T) {
	a := Vec4{1, 3, 3, 3}
	b := Vec4{1, 3, 3, 3}
	if !a.Add(b).Equals(Vec4{2, 6, 6, 6}) {
		t.Fail()
	}

	if !a.AddScalar(2).Equals(Vec4{3, 5, 5, 5}) {
		t.Fail()
	}
}

func TestVec4Sub(t *testing.T) {
	a := Vec4{1, 3, 3, 3}
	b := Vec4{1, 3, 3, 3}
	if !a.Sub(b).Equals(Vec4{0, 0, 0, 0}) {
		t.Fail()
	}

	if !a.SubScalar(2).Equals(Vec4{-1, 1, 1, 1}) {
		t.Fail()
	}
}

func TestVec4Mul(t *testing.T) {
	a := Vec4{1, 3, 3, 3}
	b := Vec4{1, 3, 3, 3}
	if !a.Mul(b).Equals(Vec4{1, 9, 9, 9}) {
		t.Fail()
	}

	if !a.MulScalar(2).Equals(Vec4{2, 6, 6, 6}) {
		t.Fail()
	}
}

func TestVec4Div(t *testing.T) {
	a := Vec4{1, 3, 3, 3}
	b := Vec4{1, 3, 3, 3}
	if !a.Div(b).Equals(Vec4{1, 1, 1, 1}) {
		t.Fail()
	}

	if !a.DivScalar(2).Equals(Vec4{0.5, 3.0 / 2.0, 3.0 / 2.0, 3.0 / 2.0}) {
		t.Fail()
	}
}
