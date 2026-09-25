import argparse
import csv
import zipfile
from pathlib import Path
from collections import defaultdict


def parse_gtfs_time(s):
    # GTFS allows times above 24:00
    h, m, sec = map(int, s.split(":"))
    return h * 3600 + m * 60 + sec


def open_gtfs_file(gtfs_path, filename):
    gtfs_path = Path(gtfs_path)

    if gtfs_path.is_dir():
        return (gtfs_path / filename).open("r", encoding="utf-8-sig", newline="")

    if zipfile.is_zipfile(gtfs_path):
        zf = zipfile.ZipFile(gtfs_path)
        return zf.open(filename, "r")

    raise ValueError(f"Expected GTFS zip or directory: {gtfs_path}")


def read_stop_ids(gtfs_path):
    stop_to_id = {}

    with open_gtfs_file(gtfs_path, "stops.txt") as f:
        if isinstance(f.read(0), bytes):
            f = (line.decode("utf-8-sig") for line in f)

        reader = csv.DictReader(f)
        for row in reader:
            stop_id = row["stop_id"]
            if stop_id not in stop_to_id:
                stop_to_id[stop_id] = len(stop_to_id)

    return stop_to_id


def convert(gtfs_path, output_path):
    gtfs_path = Path(gtfs_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading stops from {gtfs_path}")
    stop_to_id = read_stop_ids(gtfs_path)
    print(f"stops={len(stop_to_id)}")

    trip_rows = defaultdict(list)

    print("Reading stop_times.txt")
    with open_gtfs_file(gtfs_path, "stop_times.txt") as f:
        if isinstance(f.read(0), bytes):
            f = (line.decode("utf-8-sig") for line in f)

        reader = csv.DictReader(f)

        required = {
            "trip_id",
            "arrival_time",
            "departure_time",
            "stop_id",
            "stop_sequence",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"stop_times.txt missing columns: {missing}")

        rows_seen = 0
        rows_skipped_unknown_stop = 0
        rows_skipped_bad_time = 0

        for row in reader:
            rows_seen += 1

            stop_id = row["stop_id"]
            if stop_id not in stop_to_id:
                rows_skipped_unknown_stop += 1
                continue

            try:
                arr = parse_gtfs_time(row["arrival_time"])
                dep = parse_gtfs_time(row["departure_time"])
                seq = int(row["stop_sequence"])
            except Exception:
                rows_skipped_bad_time += 1
                continue

            trip_rows[row["trip_id"]].append((seq, stop_id, arr, dep))

            if rows_seen % 10_000_000 == 0:
                print(f"  stop_times rows read: {rows_seen}")

    print(f"stop_times_rows_seen={rows_seen}")
    print(f"rows_skipped_unknown_stop={rows_skipped_unknown_stop}")
    print(f"rows_skipped_bad_time={rows_skipped_bad_time}")
    print(f"trips_seen={len(trip_rows)}")

    print("Building minimum-weight undirected stop-pair edges")

    edges = {}
    segment_count = 0
    skipped_nonpositive = 0
    skipped_self_loops = 0

    for trip_id, rows in trip_rows.items():
        rows.sort(key=lambda x: x[0])

        for a, b in zip(rows, rows[1:]):
            _, stop_a, _arr_a, dep_a = a
            _, stop_b, arr_b, _dep_b = b

            u = stop_to_id[stop_a]
            v = stop_to_id[stop_b]

            if u == v:
                skipped_self_loops += 1
                continue

            w = arr_b - dep_a
            if w <= 0:
                skipped_nonpositive += 1
                continue

            x, y = (u, v) if u < v else (v, u)
            old = edges.get((x, y))
            if old is None or w < old:
                edges[(x, y)] = w

            segment_count += 1

    print(f"trip_segments_seen={segment_count}")
    print(f"distinct_undirected_edges={len(edges)}")
    print(f"skipped_self_loops={skipped_self_loops}")
    print(f"skipped_nonpositive_segments={skipped_nonpositive}")

    print(f"Writing weighted edge list to {output_path}")
    with output_path.open("w") as out:
        out.write("# Weighted GTFS stop graph\n")
        out.write("# vertex = stop_id mapped to dense integer id\n")
        out.write("# edge = consecutive stops in at least one trip\n")
        out.write("# weight = minimum scheduled in-vehicle travel time in seconds\n")
        for (u, v), w in sorted(edges.items()):
            out.write(f"{u} {v} {w}\n")

    mapping_path = output_path.with_suffix(output_path.suffix + ".stop_id_map")
    print(f"Writing stop id map to {mapping_path}")
    with mapping_path.open("w") as out:
        out.write("dense_id stop_id\n")
        for stop_id, dense_id in sorted(stop_to_id.items(), key=lambda x: x[1]):
            out.write(f"{dense_id} {stop_id}\n")

    print("done")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gtfs_zip_or_dir")
    parser.add_argument("output_edges")
    args = parser.parse_args()

    convert(args.gtfs_zip_or_dir, args.output_edges)


if __name__ == "__main__":
    main()
