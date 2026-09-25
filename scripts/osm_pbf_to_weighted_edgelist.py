import argparse
import math
from pathlib import Path

import osmium


EXCLUDED_HIGHWAYS = {
    "footway",
    "cycleway",
    "path",
    "steps",
    "pedestrian",
    "bridleway",
    "construction",
    "proposed",
    "raceway",
    "corridor",
    "elevator",
    "platform",
}


def haversine_meters(lat1, lon1, lat2, lon2):
    r = 6371008.8
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    )
    return 2.0 * r * math.asin(math.sqrt(a))


def is_road_way(way):
    highway = way.tags.get("highway")
    if highway is None:
        return False
    if highway in EXCLUDED_HIGHWAYS:
        return False
    if way.tags.get("area") == "yes":
        return False
    return True


class RoadEdgeListWriter(osmium.SimpleHandler):
    def __init__(self, out_path, directed):
        super().__init__()
        self.out_path = Path(out_path)
        self.directed = directed
        self.id_map = {}
        self.next_id = 0
        self.edges_written = 0
        self.ways_seen = 0
        self.road_ways_seen = 0
        self.bad_locations = 0
        self.zero_length_edges = 0

        self.out_path.parent.mkdir(parents=True, exist_ok=True)
        self.out = self.out_path.open("w")

        self.out.write("# weighted edge list extracted from OSM PBF\n")
        self.out.write("# format: src dst weight_in_integer_meters\n")

    def dense_id(self, osm_node_id):
        x = self.id_map.get(osm_node_id)
        if x is None:
            x = self.next_id
            self.id_map[osm_node_id] = x
            self.next_id += 1
        return x

    def way(self, way):
        self.ways_seen += 1

        if not is_road_way(way):
            return

        self.road_ways_seen += 1

        nodes = list(way.nodes)
        if len(nodes) < 2:
            return

        oneway = way.tags.get("oneway")
        is_oneway = oneway in {"yes", "true", "1"}

        for a, b in zip(nodes, nodes[1:]):
            try:
                if not a.location.valid() or not b.location.valid():
                    self.bad_locations += 1
                    continue

                dist = haversine_meters(
                    a.location.lat,
                    a.location.lon,
                    b.location.lat,
                    b.location.lon,
                )
            except osmium.InvalidLocationError:
                self.bad_locations += 1
                continue

            w = int(round(dist))
            if w <= 0:
                self.zero_length_edges += 1
                w = 1

            u = self.dense_id(a.ref)
            v = self.dense_id(b.ref)

            if u == v:
                continue

            self.out.write(f"{u} {v} {w}\n")
            self.edges_written += 1

            if self.directed and not is_oneway:
                self.out.write(f"{v} {u} {w}\n")
                self.edges_written += 1

    def close(self):
        self.out.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pbf")
    parser.add_argument("output_edges")
    parser.add_argument(
        "--directed",
        action="store_true",
        help="Emit directed edges. Non-oneway roads get both directions; oneway roads get one direction.",
    )
    args = parser.parse_args()

    handler = RoadEdgeListWriter(args.output_edges, args.directed)

    handler.apply_file(args.input_pbf, locations=True)

    handler.close()

    print(f"input_pbf={args.input_pbf}")
    print(f"output_edges={args.output_edges}")
    print(f"dense_vertices_seen={handler.next_id}")
    print(f"ways_seen={handler.ways_seen}")
    print(f"road_ways_seen={handler.road_ways_seen}")
    print(f"edge_entries_written={handler.edges_written}")
    print(f"bad_locations={handler.bad_locations}")
    print(f"zero_length_edges_clamped_to_one={handler.zero_length_edges}")


if __name__ == "__main__":
    main()
