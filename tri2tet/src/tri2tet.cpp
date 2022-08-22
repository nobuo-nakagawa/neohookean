#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <algorithm> 

#include "bary.h" // from https://www.cdsimpson.net/2014/10/barycentric-coordinates.html
#pragma warning(disable:4996) // for sscanf

// usage
// tri2tet.exe in.mesh out.h out.txt

#define OUTPUT_H (0) // output .h file for built-in

struct Vec3 {
  double x;
  double y;
  double z;
};

struct Tetrahedron {
  int id0;
  int id1;
  int id2;
  int id3;
};

struct Triangle {
  int id0;
  int id1;
  int id2;
};

struct VisVerts {
  int  tet_id;
  Vec3 bary;
};

struct Edge {
  int id0;
  int id1;
};

struct Context {
  std::vector<Vec3>        tet_vertices;
  std::vector<Tetrahedron> tetrahedra;
  std::vector<Vec3>        tri_vertices;
  std::vector<VisVerts>    vis_verts;
  std::vector<Triangle>    triangles;
  std::string              in_vol_tet; // .mesh
  std::string              in_srf_tri; // .obj
  std::string              out_txt;    // .txt
};

void usage() {
  printf("usage:\n"
         "tri2tet input_tet.{mesh,obj} input_tri.obj output.txt\n");
}

void process_args(Context* ctx, int argc, char* argv[]) {
  if (argc == 4) {
    ctx->in_vol_tet = std::string(argv[1]); // tet.mesh
    ctx->in_srf_tri = std::string(argv[2]); // tri.obj
    ctx->out_txt    = std::string(argv[3]); // out.txt
  }
}

void parse_tet_mesh(Context* ctx) {
  std::string s;
  auto& tet_vertices = ctx->tet_vertices;
  auto& tetrahedra   = ctx->tetrahedra;
  std::ifstream fs_vol_tet(ctx->in_vol_tet);
  if (!fs_vol_tet) {
    exit(0); // can't open
  }
  while(std::getline(fs_vol_tet, s)) { // Vertices
    if (s.compare("Vertices") == 0) {
      std::getline(fs_vol_tet, s);
      int num;
      if (sscanf(s.c_str(), "%d", &num) == 1) { // 320
        tet_vertices.resize(num);
      }
      for(int i = 0; i < num; i++) {
        std::getline(fs_vol_tet, s);
        double dummy; // 0! (I don't know)
        int ret = sscanf(s.c_str(), "%lf %lf %lf %lf", &tet_vertices[i].x, &tet_vertices[i].y, &tet_vertices[i].z, &dummy);
      }
    }
    if (s.compare("Tetrahedra") == 0) { // Tetrahedra
      std::getline(fs_vol_tet, s);
      int num;
      if (sscanf(s.c_str(), "%d", &num) == 1) { // 1144
        tetrahedra.resize(num);
      }
      for(int i = 0; i < num; i++) {
        std::getline(fs_vol_tet, s);
        int dummy; // 0! (I don't know)
        if (sscanf(s.c_str(), "%d %d %d %d %d", &tetrahedra[i].id0, &tetrahedra[i].id1, &tetrahedra[i].id2, &tetrahedra[i].id3, &dummy) == 5) {
          tetrahedra[i].id0 -= 1; // 1 origin to 0 origin
          tetrahedra[i].id1 -= 1; // 1 origin to 0 origin
          tetrahedra[i].id2 -= 1; // 1 origin to 0 origin
          tetrahedra[i].id3 -= 1; // 1 origin to 0 origin
        }
      }
    }
  }
}

void parse_tet_obj(Context* ctx) {
  std::string s;
  auto& tet_vertices = ctx->tet_vertices;
  auto& tetrahedra   = ctx->tetrahedra;
  std::ifstream fs_vol_tet(ctx->in_vol_tet);
  if (!fs_vol_tet) {
    exit(0); // can't open
  }
  while(std::getline(fs_vol_tet, s)) { // v 0.0 1.0 1.0
    Vec3 v;
    int ret = sscanf(s.c_str(), "v %lf %lf %lf", &v.x, &v.y, &v.z);
    if (ret == 3) {
      ctx->tet_vertices.push_back(v);
    }
    Tetrahedron tet;
    int f0, f1, f2, f3, f4, f5, f6, f7;
    ret = sscanf(s.c_str(), "f %d//%d %d//%d %d//%d %d//%d", &f0, &f1, &f2, &f3, &f4, &f5, &f6, &f7); // f 0//1 1//1 2//2 3//3
    if (ret == 8) {
      tet.id0 = f0 - 1;
      tet.id1 = f2 - 1;
      tet.id2 = f4 - 1;
      tet.id3 = f6 - 1;
      ctx->tetrahedra.push_back(tet);
    }
  }
}

void read_tet(Context* ctx) {
  if (ctx->in_vol_tet.ends_with(".mesh")) {
    parse_tet_mesh(ctx); // tet.mesh
  } else if (ctx->in_vol_tet.ends_with(".obj")) {
    parse_tet_obj(ctx);  // tet.obj
  }
}

void read_tri(Context* ctx) {
  std::string s;
  std::ifstream fs_srf_tri(ctx->in_srf_tri);
  if (!fs_srf_tri) {
    exit(0);
  }
  while(std::getline(fs_srf_tri, s)) { // v 0.0 1.0 1.0
    Vec3 v;
    int ret = sscanf(s.c_str(), "v %lf %lf %lf", &v.x, &v.y, &v.z);
    if (ret == 3) {
      ctx->tri_vertices.push_back(v);
    }
    Triangle tri;
    int f0, f1, f2, f3, f4, f5;
    ret = sscanf(s.c_str(), "f %d//%d %d//%d %d//%d", &f0, &f1, &f2, &f3, &f4, &f5); // f 0//1 1//1 2//2
    if (ret == 6) {
      tri.id0 = f0 - 1;
      tri.id1 = f2 - 1;
      tri.id2 = f4 - 1;
      ctx->triangles.push_back(tri);
    }
    else {
      ret = sscanf(s.c_str(), "f %d %d %d", &tri.id0, &tri.id1, &tri.id2); // f 0 1 2
      if (ret == 3) {
        tri.id0 -= 1; // 1 origin to 0 origin
        tri.id1 -= 1;
        tri.id2 -= 1;
        ctx->triangles.push_back(tri);
      }
    }
  }
}

void calc_weight(Context* ctx) {
  double ans[4];
  double p[3];
  double a[3];
  double b[3];
  double c[3];
  double d[3];
  for(auto& v : ctx->tri_vertices) {
    p[0] = v.x;  p[1] = v.y; p[2] = v.z;
    int idx = 0;
    for(auto& t : ctx->tetrahedra) {
      Vec3 va = ctx->tet_vertices[t.id0];
      Vec3 vb = ctx->tet_vertices[t.id1];
      Vec3 vc = ctx->tet_vertices[t.id2];
      Vec3 vd = ctx->tet_vertices[t.id3];
      a[0] = va.x;  a[1] = va.y; a[2] = va.z;
      b[0] = vb.x;  b[1] = vb.y; b[2] = vb.z;
      c[0] = vc.x;  c[1] = vc.y; c[2] = vc.z;
      d[0] = vd.x;  d[1] = vd.y; d[2] = vd.z;
      bary_tet(&ans[0], p, a, b, c, d);
      if ((ans[0] > -FLT_EPSILON) && (ans[1] > -FLT_EPSILON) && (ans[2] > -FLT_EPSILON) && (ans[3] > -FLT_EPSILON)){ // inside tetrahedron
        VisVerts vis_vert;
        vis_vert.tet_id = idx;
        vis_vert.bary.x = ans[0];
        vis_vert.bary.y = ans[1];
        vis_vert.bary.z = ans[2];
        ctx->vis_verts.push_back(vis_vert);
        break;
      }
      idx++;
    }
  }
}

void output_h(Context* ctx){
  FILE* fp = fopen(ctx->out_txt.c_str(), "w");
  if (!fp) { exit(0); } // can't open
  size_t vtx_num = ctx->tet_vertices.size();
  fprintf(fp, "std::array<Float, %d> dragonTetVerts = {\n", (int)vtx_num * 3); // tet_verts
  for(auto& v : ctx->tet_vertices) {
    fprintf(fp, "  %lf, %lf, %lf,\n", v.x, v.y, v.z);
  }
  fprintf(fp, "};\n");
  size_t tet_num = ctx->tetrahedra.size();
  fprintf(fp, "std::array<int, %d> dragonTetIds = {\n", (int)tet_num * 4); // tet_ids
  for(auto& t : ctx->tetrahedra) {
    fprintf(fp, "  %d, %d, %d, %d,\n", t.id0, t.id1, t.id2, t.id3);
  }
  fprintf(fp, "};\n");
  fprintf(fp, "std::array<int, %d> dragonTetEdgeIds = {\n", (int)tet_num * 6 * 2); // edge_id
  for(auto& t : ctx->tetrahedra) {
    fprintf(fp, "  %d, %d,\n", t.id0, t.id1);
    fprintf(fp, "  %d, %d,\n", t.id1, t.id2);
    fprintf(fp, "  %d, %d,\n", t.id2, t.id0);
    fprintf(fp, "  %d, %d,\n", t.id0, t.id3);
    fprintf(fp, "  %d, %d,\n", t.id1, t.id3);
    fprintf(fp, "  %d, %d,\n", t.id2, t.id3);
  }
  fprintf(fp, "};\n");
  fprintf(fp, "std::array<Float, %d> dragonAttachedVerts = {\n", (int)ctx->vis_verts.size() * 4); // vis verts
  for(auto& t : ctx->vis_verts) {
    fprintf(fp, "  %d, %lf, %lf, %lf,\n", t.tet_id, t.bary.x, t.bary.y, t.bary.z);
  }
  fprintf(fp, "};\n");
  size_t tri_num = ctx->triangles.size();
  fprintf(fp, "std::array<int, %d> dragonAttachedTriIds = {\n", (int)tri_num * 3); // vis tri id
  for(auto& t : ctx->triangles) {
    fprintf(fp, "  %d, %d, %d,\n", t.id0, t.id1, t.id2);
  }
  fprintf(fp, "};\n");
  fclose(fp);
}

void output_txt(Context* ctx) {
  FILE* fp = fopen(ctx->out_txt.c_str(), "w");
  if (!fp) { exit(0); } // can't open
  fprintf(fp, "<tetrahedron vertices %d>\n", (int)ctx->tet_vertices.size());
  for(auto& v : ctx->tet_vertices) {
    fprintf(fp, "  %lf, %lf, %lf,\n", v.x, v.y, v.z);
  }
  fprintf(fp, "<tetrahedron id %d>\n", (int)ctx->tetrahedra.size());
  for(auto& t : ctx->tetrahedra) {
    fprintf(fp, "  %d, %d, %d, %d,\n", t.id0, t.id1, t.id2, t.id3);
  }
#if 0
  fprintf(fp, "<tetrahedron edge id %d>\n", (int)ctx->tetrahedra.size() * 6 * 2);
  for(auto& t : ctx->tetrahedra) {
    fprintf(fp, "  %d, %d,\n", t.id0, t.id1);
    fprintf(fp, "  %d, %d,\n", t.id1, t.id2);
    fprintf(fp, "  %d, %d,\n", t.id2, t.id0);
    fprintf(fp, "  %d, %d,\n", t.id0, t.id3);
    fprintf(fp, "  %d, %d,\n", t.id1, t.id3);
    fprintf(fp, "  %d, %d,\n", t.id2, t.id3);
  }
#else
  std::vector<Edge> edges;
  for(auto& t : ctx->tetrahedra) {
    Edge edge;
    edge.id0 = std::min(t.id0, t.id1); // id0-id1
    edge.id1 = std::max(t.id0, t.id1); // id0 < id1 (sorted)
    edges.push_back(edge);
    edge.id0 = std::min(t.id1, t.id2); // id1-id2
    edge.id1 = std::max(t.id1, t.id2); // id1 < id2 (sorted)
    edges.push_back(edge);
    edge.id0 = std::min(t.id2, t.id0); // id2-id0
    edge.id1 = std::max(t.id2, t.id0); // id2 < id0 (sorted)
    edges.push_back(edge);
    edge.id0 = std::min(t.id0, t.id3); // id0-id3
    edge.id1 = std::max(t.id0, t.id3); // id0 < id3 (sorted)
    edges.push_back(edge);
    edge.id0 = std::min(t.id1, t.id3); // id1-id3
    edge.id1 = std::max(t.id1, t.id3); // id1 < id3 (sorted)
    edges.push_back(edge);
    edge.id0 = std::min(t.id2, t.id3); // id2-id3
    edge.id1 = std::max(t.id2, t.id3); // id2 < id3 (sorted)
    edges.push_back(edge);
  }
  auto comp = [](const auto& lh, const auto& rh) {
    return (lh.id0 * 65535 + lh.id1) < (rh.id0 * 65535 + rh.id1);
  };
  auto pred = [](const auto& lh, const auto& rh) {
    return (lh.id0 == rh.id0) && (lh.id1 == rh.id1);
  };
  std::sort(edges.begin(), edges.end(), comp);
  edges.erase(std::unique(edges.begin(), edges.end(), pred), edges.end());
  fprintf(fp, "<tetrahedron edge id %d>\n", (int)edges.size() * 2);
  for(auto& e : edges) {
    fprintf(fp, "  %d, %d,\n", e.id0, e.id1); // unique
  }
#endif
  fprintf(fp, "<attached triangle vertices %d>\n", (int)ctx->vis_verts.size());
  for(auto& t : ctx->vis_verts) {
    fprintf(fp, "  %d, %lf, %lf, %lf,\n", t.tet_id, t.bary.x, t.bary.y, t.bary.z);
  }
  fprintf(fp, "<attached triangle id %d>\n", (int)ctx->triangles.size());
  for(auto& t : ctx->triangles) {
    fprintf(fp, "  %d, %d, %d,\n", t.id0, t.id1, t.id2);
  }
  fclose(fp);
}

int main(int argc, char* argv[]) {
  Context context;
  if (argc != 4) {
    usage();
  } else {
    process_args(&context, argc, argv);
    read_tet(&context);
    read_tri(&context);
    calc_weight(&context);
#if OUTPUT_H
    output_h(&context);
#else
    output_txt(&context);
#endif
  }
}
