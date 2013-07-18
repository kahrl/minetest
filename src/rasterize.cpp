/*
Minetest
Copyright (C) 2010-2013 celeron55, Perttu Ahola <celeron55@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "rasterize.h"
#include "irrlichttypes_extrabloated.h"
#include "debug.h"
#include "log.h"
#include "util/numeric.h"
#include "util/timetaker.h"
#include <algorithm> // std::reverse

#define PP2(x) "("<<(x).X<<","<<(x).Y<<")"
#define PP(x) "("<<(x).X<<","<<(x).Y<<","<<(x).Z<<")"

bool make_convex_polygon(ConvexPolygon &pol, const std::vector<v2f> &vertices)
{
	u32 k;

	pol.vl.clear();
	pol.vr.clear();

	if(vertices.size() <= 2)
		return false;

  	// find index i of vertex such that its y coordinate is minimal
	u32 i = 0;
	for(k = 1; k < vertices.size(); ++k){
		if(vertices[k].Y < vertices[i].Y)
			i = k;
	}

	// we determined the topmost vertex; rotate the rest of the vertex list
	v2f topmost = vertices[i];
	std::vector<v2f> rest;
	rest.insert(rest.end(), vertices.begin()+i+1, vertices.end());
	rest.insert(rest.end(), vertices.begin(), vertices.begin()+i);

	// in case vertices were not given in CCW order, return false
	v2f anglevec1 = rest.front() - topmost;
	v2f anglevec2 = rest.back() - topmost;
	// use cross product to determine sign of angle between
	// anglevec1 and anglevec2
	if(anglevec1.X*anglevec2.Y - anglevec1.Y*anglevec2.X > 0){
		dstream<<"clockwise vertices given"<<std::endl;
		return false;
	}

	// find smallest index j such that rest[j].y >= rest[j+1].y
	// or j = rest.size()-1 if such an index does not exist
	u32 j = 0;
	while(j+1 < rest.size() && rest[j].Y < rest[j+1].Y)
		++j;

	// split rest into vl (before j) and vr (after j) such that vl and vr
	// begin with topmost == vertices[i] and end with bottommost == rest[j]
	pol.vl.push_back(topmost);
	pol.vl.insert(pol.vl.end(), rest.begin(), rest.begin()+j+1);
	pol.vr.insert(pol.vr.end(), rest.begin()+j, rest.end());
	pol.vr.push_back(topmost);
	std::reverse(pol.vr.begin(), pol.vr.end());

	// remove horizontal line segments from the top of the polygon
	k = 0;
	while(pol.vl.size() >= k+2 && pol.vl[k].Y == pol.vl[k+1].Y)
		++k;
	pol.vl.erase(pol.vl.begin(), pol.vl.begin() + k);
	k = 0;
	while(pol.vr.size() >= k+2 && pol.vr[k].Y == pol.vr[k+1].Y)
		++k;
	pol.vr.erase(pol.vr.begin(), pol.vr.begin() + k);

	// remove horizontal line segments from the bottom of the polygon
	while(pol.vl.size() >= 2 &&
			pol.vl[pol.vl.size()-2].Y == pol.vl[pol.vl.size()-1].Y){
		pol.vl.pop_back();
	}
	while(pol.vr.size() >= 2 &&
			pol.vr[pol.vr.size()-2].Y == pol.vr[pol.vr.size()-1].Y){
		pol.vr.pop_back();
	}

	return true;
}

void dump_convex_polygon(const ConvexPolygon &pol)
{
	for(u32 i = 0; i < pol.vl.size(); ++i){
		dstream<<"vl["<<i<<"] = ("<<pol.vl[i].X<<","<<pol.vl[i].Y<<")"<<std::endl;
	}
	for(u32 i = 0; i < pol.vr.size(); ++i){
		dstream<<"vr["<<i<<"] = ("<<pol.vr[i].X<<","<<pol.vr[i].Y<<")"<<std::endl;
	}
}

static void rasterize_interpolate_x(v2f p1, v2f p2, f32& slope, f32 &intercept)
{
	slope = (p2.Y - p1.Y) / (p2.X - p1.X);
	intercept = p1.Y - p1.X * slope;
}

static void rasterize_interpolate_y(v2f p1, v2f p2, f32& slope, f32 &intercept)
{
	slope = (p2.X - p1.X) / (p2.Y - p1.Y);
	intercept = p1.X - p1.Y * slope;
}

static void draw_convex_polygon(video::IImage *img_src, video::IImage *img_dst,
		const video::SColor &color, const ConvexPolygon &pol,
		const core::matrix4& T)
{
	if(pol.vl.size() <= 1 || pol.vr.size() <= 1)
		return;

	core::dimension2d<u32> dstdim = img_dst->getDimension();
	u32 w = dstdim.Width;
	u32 h = dstdim.Height;

	//void *data = img_dst->lock();
	//if(data == NULL)
	//	return;

	// Clipping rectangle (clip_ymin is updated dynamically)
	s32 clip_xmin = 0;
	s32 clip_xmax = w - 1;
	s32 clip_ymin = 0;
	s32 clip_ymax = h - 1;

	// Split polygon into trapezoids and draw them, starting at the top
	v2f pl = pol.vl[0];
	v2f pr = pol.vr[0];
	u32 jl = 1;
	u32 jr = 1;
	v2f pl2, pr2;
	int jl2, jr2;

	while(jl < pol.vl.size() && jr < pol.vr.size()){
		//dstream<<"jl = "<<jl<<", jr="<<jr<<std::endl;

		// Find the next trapezoid
		if(pol.vl[jl].Y <= pol.vr[jr].Y){
			pl2 = pol.vl[jl];
			jl2 = jl + 1;
			f32 slope, intercept;
			rasterize_interpolate_y(pol.vr[jr-1], pol.vr[jr],
					slope, intercept);
			pr2 = v2f(slope * pl2.Y + intercept, pl2.Y);
			jr2 = jr;
		}
		else{
			pr2 = pol.vr[jr];
			jr2 = jr + 1;
			f32 slope, intercept;
			rasterize_interpolate_y(pol.vl[jl-1], pol.vl[jl],
					slope, intercept);
			pl2 = v2f(slope * pr2.Y + intercept, pr2.Y);
			jl2 = jl;
		}

		// Note that pl.Y == pr.Y and pl2.Y == pr2.Y
		// Draw trapezoid bounded by the corners pl, pl2, pr2, pr
		#if 0
		dstream<<"Drawing trapezoid:"<<std::endl;
		dstream<<"pl  = "<<PP2(pl)<<std::endl;
		dstream<<"pl2 = "<<PP2(pl2)<<std::endl;
		dstream<<"pr  = "<<PP2(pr)<<std::endl;
		dstream<<"pr2 = "<<PP2(pr2)<<std::endl;
		#endif
		s32 ymin = MYMAX(floor(pl.Y), clip_ymin);
		s32 ymax = MYMIN(ceil(pl2.Y), clip_ymax);
		f32 xmin_slope, xmin_intercept;
		f32 xmax_slope, xmax_intercept;
		rasterize_interpolate_y(pl, pl2, xmin_slope, xmin_intercept);
		rasterize_interpolate_y(pr, pr2, xmax_slope, xmax_intercept);
		for(s32 y = ymin; y <= ymax; ++y){
			s32 xmin = MYMAX(floor(xmin_slope*y + xmin_intercept),
					clip_xmin);
			s32 xmax = MYMIN(ceil(xmax_slope*y + xmax_intercept),
					clip_xmax);
			#if 0
			dstream<<"Drawing scanline "<<y<<": "<<xmin<<" .. "<<xmax<<std::endl;
			#endif
			for(s32 x = xmin; x <= xmax; ++x){
				// Get source texture coordinates
				s32 tx = T[0]*x + T[4]*y + T[12] + 0.5;
				s32 ty = T[1]*x + T[5]*y + T[13] + 0.5;
				// Get Z-buffer value
				//s32 z = T[2]*x + T[6] + T[14];
				img_dst->setPixel(x, y, img_src->getPixel(tx, ty));
			}
		}

		// Update clipping rect to avoid drawing to a scanline twice
		clip_ymin = ymax + 1;

		// Update loop variables
		pl = pl2;
		pr = pr2;
		jl = jl2;
		jr = jr2;
	}
}

static void draw_3d_quad(video::IImage *img_src, video::IImage *img_dst,
		const video::SColor &color, const core::matrix4 A,
		v3f w00, v3f w01, v3f w11, v3f w10)
{
	//u32 *data = (u32*) img_dst->lock();
	//if(data == NULL)
	//	return;
	
	v3f p00, p01, p11, p10;
	A.transformVect(p00, w00);
	A.transformVect(p01, w01);
	A.transformVect(p11, w11);
	A.transformVect(p10, w10);
	#if 0
	dstream<<"p00="<<PP(p00)<<std::endl;
	dstream<<"p01="<<PP(p01)<<std::endl;
	dstream<<"p11="<<PP(p11)<<std::endl;
	dstream<<"p10="<<PP(p10)<<std::endl;
	#endif

	f32 wsrc = img_src->getDimension().Width - 1;
	f32 hsrc = img_src->getDimension().Height - 1;

	// generate a matrix T such that (px,py,0,1)*A = (tx,ty,z,1)
	// where (px,py) = projection space coordinates,
	//       (tx,ty) = texture coordinates (absolute in img_src),
	//       z       = z buffer value
	//
	// in particular:
	// (p00.X, p00.Y, 0, 1)*T = (0,       0, p00.Z, 1)
	// (p10.X, p10.Y, 0, 1)*T = (wsrc,    0, p10.Z, 1)
	// (p01.X, p01.Y, 0, 1)*T = (0,    hsrc, p01.Z, 1)

	core::matrix4 T1;
	T1.setTranslation(v3f(-p00.X, -p00.Y, 0));

	v3f q10 = p10 - p00;
	v3f q01 = p01 - p00;
	f32 t2_invdet = q10.X * q01.Y - q10.Y * q01.X;
	if(fabs(t2_invdet) < 1e-10){
		verbosestream<<"draw_3d_quad: unable to compute texture matrix,"
			<<"t2_invdet is zero"<<std::endl;
		return;
	}
	f32 t2_det = 1.0 / t2_invdet;
	core::matrix4 T2;
	T2[0] = +q01.Y * t2_det;
	T2[1] = -q10.Y * t2_det;
	T2[4] = -q01.X * t2_det;
	T2[5] = +q10.X * t2_det;
	T2[10] = t2_det;

	core::matrix4 T3;
	T3[0] = wsrc;
	T3[2] = p10.Z - p00.Z;
	T3[5] = hsrc;
	T3[6] = p01.Z - p00.Z;
	T3[14] = p00.Z;

	// note: irrlicht's matrix multiplication order is the reverse
	// of the common multiplication order used in mathematics
	core::matrix4 T = T3 * T2 * T1;

	#if 0
	dstream<<"wsrc = "<<wsrc<<", hsrc = "<<hsrc<<std::endl;
	dstream<<"p00 = "<<PP(p00)<<std::endl;
	dstream<<"p10 = "<<PP(p10)<<std::endl;
	dstream<<"p01 = "<<PP(p01)<<std::endl;
	v3f tmp00, tmp10, tmp01;
	T.transformVect(tmp00, v3f(p00.X, p00.Y, 0));
	T.transformVect(tmp10, v3f(p10.X, p10.Y, 0));
	T.transformVect(tmp01, v3f(p01.X, p01.Y, 0));
	dstream<<"((p00.X, p00.Y, 0), 1)*T = ("<<PP(tmp00)<<",1)"<<std::endl;
	dstream<<"((p10.X, p10.Y, 0), 1)*T = ("<<PP(tmp10)<<",1)"<<std::endl;
	dstream<<"((p01.X, p01.Y, 0), 1)*T = ("<<PP(tmp01)<<",1)"<<std::endl;
	#endif

	// generate a convex polygon that describes the area to be rasterized
	ConvexPolygon pol;
	std::vector<v2f> vertices;
	vertices.push_back(v2f(p00.X, p00.Y));
	vertices.push_back(v2f(p01.X, p01.Y));
	vertices.push_back(v2f(p11.X, p11.Y));
	vertices.push_back(v2f(p10.X, p10.Y));
	make_convex_polygon(pol, vertices);
	//dump_convex_polygon(pol);
	draw_convex_polygon(img_src, img_dst, color, pol, T);
}

void inventorycube(video::IImage *img_top, video::IImage *img_left,
		video::IImage *img_right, video::IImage *img_dst)
{
	TimeTaker timetaker("inventorycube", NULL, PRECISION_NANO);

	if(img_top == NULL || img_left == NULL || img_right == NULL
			|| img_dst == NULL)
		return;

	assert(img_dst->getColorFormat() == video::ECF_A8R8G8B8);

	core::dimension2d<u32> dstdim = img_dst->getDimension();
	u32 w = dstdim.Width;
	u32 h = dstdim.Height;

	void *data = img_dst->lock();
	if(data == NULL)
		return;

	// Clear destination image
	memset(data, 0, img_dst->getImageDataSizeInBytes());

	// Calculate projection-to-view transformation
	// (orthogonal, left-handed)

	f32 pw = 1.65;  // width of view volume
	f32 ph = 1.65;  // height of view volume
	v2f pscale(pw / w, -ph / h);
	v2f ptrans(-pw / 2, ph / 2);

	// Calculate view-to-world transformation (left-handed)

	v3f position(0, 1.0, -1.5);  // camera position
	position.rotateXZBy(45);     // ^
	v3f target(0, 0, 0);         // camera look-at point
	v3f upVector(0, 1, 0);       // camera up vector

	v3f zaxis = target - position;
	zaxis.normalize();

	v3f xaxis = upVector.crossProduct(zaxis);
	xaxis.normalize();

	v3f yaxis = zaxis.crossProduct(xaxis);

	// Calculate projection-to-world transformation

	core::matrix4 pt;  // (part of a) matrix

	pt[0]  = pscale.X*xaxis.X;
	pt[1]  = pscale.X*xaxis.Y;
	pt[2]  = pscale.X*xaxis.Z;

	pt[4]  = pscale.Y*yaxis.X;
	pt[5]  = pscale.Y*yaxis.Y;
	pt[6]  = pscale.Y*yaxis.Z;

	pt[8]  = zaxis.X;
	pt[9]  = zaxis.Y;
	pt[10] = zaxis.Z;

	pt[12] = ptrans.X * xaxis.X + ptrans.Y * yaxis.X + position.X;
	pt[13] = ptrans.X * xaxis.Y + ptrans.Y * yaxis.Y + position.Y;
	pt[14] = ptrans.X * xaxis.Z + ptrans.Y * yaxis.Z + position.Z;

	core::matrix4 A;
	pt.getInverse(A);

	// "Lighting".
	video::SColor color_top(255, 255, 255, 255);
	video::SColor color_left(255, 255, 255, 255);
	video::SColor color_right(255, 255, 255, 255);

	//dstream<<"draw_3d_quad(img_top)"<<std::endl;
	draw_3d_quad(img_top, img_dst, color_top, A,
			v3f(-0.5,+0.5,+0.5), v3f(-0.5,+0.5,-0.5),
			v3f(+0.5,+0.5,-0.5), v3f(+0.5,+0.5,+0.5));
	//dstream<<"draw_3d_quad(img_left)"<<std::endl;
	draw_3d_quad(img_left, img_dst, color_left, A,
			v3f(-0.5,+0.5,-0.5), v3f(-0.5,-0.5,-0.5),
			v3f(+0.5,-0.5,-0.5), v3f(+0.5,+0.5,-0.5));
	//dstream<<"draw_3d_quad(img_right)"<<std::endl;
	draw_3d_quad(img_right, img_dst, color_right, A,
			v3f(+0.5,+0.5,-0.5), v3f(+0.5,-0.5,-0.5),
			v3f(+0.5,-0.5,+0.5), v3f(+0.5,+0.5,+0.5));

	timetaker.stop(true);
	infostream<<"inventorycube took "<<timetaker.getTimerTime()<<"ns"<<std::endl;
}

