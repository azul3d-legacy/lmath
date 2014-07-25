// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestVec2Add(t *testing.T) {
	a := Vec2{1, 3}
	b := Vec2{1, 3}
	if !a.Add(b).Equals(Vec2{2, 6}) {
		t.Fail()
	}

	if !a.AddScalar(2).Equals(Vec2{3, 5}) {
		t.Fail()
	}
}

func TestVec2Sub(t *testing.T) {
	a := Vec2{1, 3}
	b := Vec2{1, 3}
	if !a.Sub(b).Equals(Vec2{0, 0}) {
		t.Fail()
	}

	if !a.SubScalar(2).Equals(Vec2{-1, 1}) {
		t.Fail()
	}
}

func TestVec2Mul(t *testing.T) {
	a := Vec2{1, 3}
	b := Vec2{1, 3}
	if !a.Mul(b).Equals(Vec2{1, 9}) {
		t.Fail()
	}

	if !a.MulScalar(2).Equals(Vec2{2, 6}) {
		t.Fail()
	}
}

func TestVec2Div(t *testing.T) {
	a := Vec2{1, 3}
	b := Vec2{1, 3}
	if !a.Div(b).Equals(Vec2{1, 1}) {
		t.Fail()
	}

	if !a.DivScalar(2).Equals(Vec2{0.5, 3.0 / 2.0}) {
		t.Fail()
	}
}
