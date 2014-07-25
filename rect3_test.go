// Copyright 2014 The Azul3D Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package lmath

import (
	"testing"
)

func TestRect3SqDistToPoint(t *testing.T) {
	r := Rect3{
		Min: Vec3{-1, -1, -1},
		Max: Vec3{1, 1, 1},
	}
	p := Vec3{-2, 2, -2}
	wantDist2 := r.Closest(p).Sub(p).LengthSq()
	gotDist2 := r.SqDistToPoint(p)
	if !Equal(wantDist2, gotDist2) {
		t.Log("want", wantDist2)
		t.Log("got", gotDist2)
		t.Fail()
	}
}

func BenchmarkRect3SqDistToPointSlow(b *testing.B) {
	r := Rect3{
		Min: Vec3{-1, -1, -1},
		Max: Vec3{1, 1, 1},
	}
	p := Vec3{-2, 2, -2}
	for n := 0; n < b.N; n++ {
		r.Closest(p).Sub(p).LengthSq()
	}
}

func BenchmarkRect3SqDistToPointFast(b *testing.B) {
	r := Rect3{
		Min: Vec3{-1, -1, -1},
		Max: Vec3{1, 1, 1},
	}
	p := Vec3{-2, 2, -2}
	for n := 0; n < b.N; n++ {
		r.SqDistToPoint(p)
	}
}
