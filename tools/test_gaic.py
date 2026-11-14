import math
import os
import numpy as np
from openalea.plantgl.all import Scene, TriangleSet


def find_triangles_in_radius(scene_or_file, target_point, radius, xy_only=True):

    scene = Scene(scene_or_file) if isinstance(scene_or_file, str) else scene_or_file
    target = np.array(target_point[:2] if xy_only else target_point)
    results = []
    total_triangle_area = 0.0

    for shape_idx, shape in enumerate(scene):
        geom = shape.geometry
        if not isinstance(geom, TriangleSet):
            continue

        # Extract triangles points and indices
        pts = np.array([[p.x, p.y, p.z] for p in geom.pointList])
        indices = np.array([[tri[0], tri[1], tri[2]] for tri in geom.indexList])

        if len(indices) == 0:
            continue

        # Get triangle vertex
        v0 = pts[indices[:, 0]]
        v1 = pts[indices[:, 1]]
        v2 = pts[indices[:, 2]]

        # Computes barycenters (3D)
        barycenters = (v0 + v1 + v2) / 3.0

        # 2D-Projection 2D (XY) ou 3D distance
        coords = barycenters[:, :2] if xy_only else barycenters

        # 3D area of triangles : 0.5 * ||(B - A) x (C - A)||
        cross_products = np.cross(v1 - v0, v2 - v0)
        areas = 0.5 * np.linalg.norm(cross_products, axis=1)

        # Euclidean distance
        distances = np.linalg.norm(coords - target, axis=1)

        # Sélection des triangles dans le rayon
        within_radius = np.where(distances <= radius)[0]

        for tri_idx in within_radius:
            area = float(areas[tri_idx])
            total_triangle_area += area
            results.append({
                'shape_id': getattr(shape, 'id', shape_idx),
                'triangle_index': int(tri_idx),
                'barycenter': barycenters[tri_idx],
                'distance': float(distances[tri_idx]),
                'area': area
            })

        # Soil surface considered for GAIp
        circle_area = math.pi * (radius ** 2)

        # GAIp
        GAIp = total_triangle_area / circle_area if circle_area > 0 else 0.0

        metrics = {
            'total_triangle_area': total_triangle_area,
            'GAIp': GAIp
        }

    return results, metrics


if __name__ == "__main__":
    folder_path = r""
    scene_file = "scene_light_plantgl_12033924.bgeom"
    scene_path = os.path.join(folder_path, scene_file)

    target = (0.125, 0.15, 0.0)
    radius = 0.075 # meters around target

    neighbour_triangles, metrics = find_triangles_in_radius(scene_path, target_point=target, radius=radius, xy_only=True)

    print(f"Nombre de triangles trouvés dans un rayon de {radius} m : {len(neighbour_triangles)}")
    print(f"Surface totale des triangles      : {metrics['total_triangle_area']:.6f} m²")
    print(f"GAIp : {metrics['GAIp']:.4f}")