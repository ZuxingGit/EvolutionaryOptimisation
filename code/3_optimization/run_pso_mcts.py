import multiprocessing
from shapely.geometry import Polygon
from shapely.geometry import LineString
import random
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import MultiPolygon
import os
import math
from datetime import datetime
import copy

# Land boundary
min_house_width = 6
max_house_width = 18

min_house_depth = 8
max_house_depth = 25

# house_width = random.randint(min_house_width, max_house_width)
# house_depth = random.randint(min_house_depth, max_house_depth)

# set manually for now
house_width = 10
house_depth = 15

boundary = Polygon(
    [(0, 0), (house_width, 0), (house_width, house_depth), (0, house_depth)]
)

# The entry position
entry = LineString([(4.5, 0), (5.5, 0)])


# Define room types (type, min_width, max_width, min_depth, max_depth)
rooms_range = [
    ("GAR", 3.5, 6, 5.5, 7),  # Garage
    ("LDR", 2, 3, 2, 3),  # Laundry Room
    ("LR", 3, 6, 3, 6),  # Living Room
    ("MBR", 3, 6, 3, 6),  # Master Bedroom
    ("BR1", 3, 4.5, 3, 4.5),  # Bedroom
    ("KIC", 2, 4, 2, 4),  # Kitchen
    ("BA", 2, 4, 2, 5),  # Bathroom
    ("DR", 3, 5, 3, 5),  # Dining Room
    # ("BR2", 3, 4.5, 3, 4.5),  #Bedroom2
    # ("HW", 1, 1.5, 2, 10),  # Hallway
    # ("BAL", 1, 2, 2, 8),  # Balcony
]


# Components size range
min_window_width = 0.5
max_window_width = 1.5

min_door_width = 0.8
max_door_width = 1.5


class Window:
    def __init__(self, width, wall, x, y):
        self.width = width
        self.wall = wall  # left, right, top, bottom
        # window starting point(x, y)
        self.x = x
        self.y = y
        self.coordinates = self.get_coordinates()
        self.type = "WINDOW"

    def get_coordinates(self):
        if self.wall == "left" or self.wall == "right":
            return LineString([(self.x, self.y), (self.x, self.y + self.width)])
        elif self.wall == "top" or self.wall == "bottom":
            return LineString([(self.x, self.y), (self.x + self.width, self.y)])

    def get_center(self):
        if self.wall == "left" or self.wall == "right":
            return (self.x, self.y + self.width / 2)
        elif self.wall == "top" or self.wall == "bottom":
            return (self.x + self.width / 2, self.y)


class Door:
    def __init__(self, width, wall, x, y):
        self.width = width
        self.wall = wall  # left, right, top, bottom
        # door starting point(x, y)
        self.x = x
        self.y = y
        self.coordinates = self.get_coordinates()
        self.type = "DOOR"

    def get_coordinates(self):
        if self.wall == "left" or self.wall == "right":
            return LineString([(self.x, self.y), (self.x, self.y + self.width)])
        elif self.wall == "top" or self.wall == "bottom":
            return LineString([(self.x, self.y), (self.x + self.width, self.y)])

    def get_center(self):
        if self.wall == "left" or self.wall == "right":
            return (self.x, self.y + self.width / 2)
        elif self.wall == "top" or self.wall == "bottom":
            return (self.x + self.width / 2, self.y)


# Room class
class Room:
    def __init__(self, name, width, depth, x, y):
        self.name = name
        self.width = width
        self.depth = depth
        self.x = x
        self.y = y
        self.windows = []
        self.doors = []

    def get_polygon(self):
        return Polygon(
            [
                (self.x, self.y),
                (self.x + self.width, self.y),
                (self.x + self.width, self.y + self.depth),
                (self.x, self.y + self.depth),
            ]
        )

    def __str__(self):
        return f"{self.name} width: {self.width}, depth: {self.depth} at ({self.x}, {self.y})"


# House class
class House:
    def __init__(self, rooms, boundary):
        self.rooms = rooms
        self.cluster = self.get_cluster()
        self.boundary = boundary

    # combine rooms into a cluster polygon
    def get_cluster(self):
        cluster = Polygon()
        for room in self.rooms:
            cluster = cluster.union(room.get_polygon())
        return cluster

    def generate_doors_windows(self):
        for room in self.rooms:
            if room.name not in ["KIC", "LR", "DR"]:
                self.generate_door(room)
            if room.name not in ["GAR"]:
                self.generate_window(room)

    def generate_window(self, room):
        """
        generate a window for the room, no window for garage
        """
        if room.name == "GAR":
            return None

        room.windows = []
        # First check if one side of the room is facing outside
        # if not, do not generate window
        directions = ["left", "right", "top", "bottom"]
        options = []
        for direction in directions:
            if self.is_facing_outside(room, direction):
                options.append(direction)
        if not options:
            return None

        wall = random.choice(options)

        # window width, make sure not bigger than the room width or depth
        window_width = round(random.uniform(min_window_width, max_window_width), 1)
        if wall == "left" or wall == "right":
            if window_width > room.depth:
                window_width = room.depth
        elif wall == "top" or wall == "bottom":
            if window_width > room.width:
                window_width = room.width

        # window starting point
        if wall == "left":
            x = room.x
            y = round(random.uniform(room.y, room.y + room.depth - window_width), 1)
        elif wall == "right":
            x = room.x + room.width
            y = round(random.uniform(room.y, room.y + room.depth - window_width), 1)
        elif wall == "top":
            x = round(random.uniform(room.x, room.x + room.width - window_width), 1)
            y = room.y + room.depth
        elif wall == "bottom":
            x = round(random.uniform(room.x, room.x + room.width - window_width), 1)
            y = room.y

        room.windows.append(Window(window_width, wall, x, y))

    def generate_door(self, room):
        """
        generate a door for the room.
        no door for kitchen, living room, dining room.
        Garage can have multiple doors.
        """
        if room.name in ["KIC", "LR", "DR"]:
            return None

        room.doors = []
        # First check if one side of the room is facing outside
        # if not, generate door for the room
        directions = ["left", "right", "top", "bottom"]
        options = []
        for direction in directions:
            if not self.is_facing_outside(
                room, direction
            ) and self.is_facing_inside_public(room, direction):
                options.append(direction)
        if not options:
            return None

        wall = random.choice(options)

        # door width, make sure not bigger than the room width or depth
        door_width = round(random.uniform(min_door_width, max_door_width), 1)
        if wall == "left" or wall == "right":
            if door_width > room.depth:
                door_width = room.depth
        elif wall == "top" or wall == "bottom":
            if door_width > room.width:
                door_width = room.width

        # door starting point
        if wall == "left":
            x = room.x
            y = round(random.uniform(room.y, room.y + room.depth - door_width), 1)
        elif wall == "right":
            x = room.x + room.width
            y = round(random.uniform(room.y, room.y + room.depth - door_width), 1)
        elif wall == "top":
            x = round(random.uniform(room.x, room.x + room.width - door_width), 1)
            y = room.y + room.depth
        elif wall == "bottom":
            x = round(random.uniform(room.x, room.x + room.width - door_width), 1)
            y = room.y

        room.doors.append(Door(door_width, wall, x, y))

    # check if a side of a room facing outside
    def is_facing_outside(self, room, side):
        """
        Check if a side of a room is facing outside, not other rooms
        """
        if side == "left":
            for other_room in self.rooms:
                if other_room.name != room.name:
                    if (
                        other_room.x < room.x
                        and other_room.y < room.y + room.depth
                        and other_room.y + other_room.depth > room.y
                    ):
                        return False

        elif side == "right":
            for other_room in self.rooms:
                if other_room.name != room.name:
                    if (
                        other_room.x + other_room.width > room.x + room.width
                        and other_room.y < room.y + room.depth
                        and other_room.y + other_room.depth > room.y
                    ):
                        return False

        elif side == "top":
            for other_room in self.rooms:
                if other_room.name != room.name:
                    if (
                        other_room.y + other_room.depth > room.y + room.depth
                        and other_room.x < room.x + room.width
                        and other_room.x + other_room.width > room.x
                    ):
                        return False

        elif side == "bottom":
            for other_room in self.rooms:
                if other_room.name != room.name:
                    if (
                        other_room.y < room.y
                        and other_room.x < room.x + room.width
                        and other_room.x + other_room.width > room.x
                    ):
                        return False

        return True

    # check if a side of a room facing inside and not other private rooms
    def is_facing_inside_public(self, room, side):
        """
        Check if a side of a room is facing inside and touching non-private areas/rooms, e.g. LR, DR, KIC.
        Or just touching nothing, But must face inside. Must used along with `not is_facing_outside`
        """
        if side == "left":
            wall = LineString([(room.x, room.y), (room.x, room.y + room.depth)])
            for other_room in self.rooms:
                if (
                    other_room.name != room.name
                    and other_room.name != "LR"
                    and other_room.name != "DR"
                    and other_room.name != "KIC"
                ):
                    if wall.intersects(other_room.get_polygon()):
                        return False

        elif side == "right":
            wall = LineString(
                [
                    (room.x + room.width, room.y),
                    (room.x + room.width, room.y + room.depth),
                ]
            )
            for other_room in self.rooms:
                if (
                    other_room.name != room.name
                    and other_room.name != "LR"
                    and other_room.name != "DR"
                    and other_room.name != "KIC"
                ):
                    if wall.intersects(other_room.get_polygon()):
                        return False

        elif side == "top":
            wall = LineString(
                [
                    (room.x, room.y + room.depth),
                    (room.x + room.width, room.y + room.depth),
                ]
            )
            for other_room in self.rooms:
                if (
                    other_room.name != room.name
                    and other_room.name != "LR"
                    and other_room.name != "DR"
                    and other_room.name != "KIC"
                ):
                    if wall.intersects(other_room.get_polygon()):
                        return False

        elif side == "bottom":
            wall = LineString([(room.x, room.y), (room.x + room.width, room.y)])
            for other_room in self.rooms:
                if (
                    other_room.name != room.name
                    and other_room.name != "LR"
                    and other_room.name != "DR"
                    and other_room.name != "KIC"
                ):
                    if wall.intersects(other_room.get_polygon()):
                        return False

        return True


# Feasibility check
def is_valid_placement(rect, cluster, boundary):
    """
    check if the rectangle is within the house land
    but unused in the current implementation
    """
    if not boundary.contains(rect):
        return False

    # check if the rectangle attaches to the cluster, neither overlap nor separate
    if cluster is None:
        return True
    elif cluster.touches(rect):
        touch_part = cluster.intersection(rect)
        # not just a point touching
        if touch_part.geom_type in [
            "Point",
            "MultiPoint",
            "LineString",
            "MultiLineString",
        ]:
            return True
        else:
            return False
    else:
        return False


# Fitness functions


## 1. dis(MBR, BR) +
# calculate the distance between the center of the master bedroom and bedroom
# normalize the distance
def dis_MBR_BR(house):
    MBR = None
    BR = None
    for room in house.rooms:
        if room.name == "MBR":
            MBR = room
        if room.name.startswith("BR"):
            BR = room
    if MBR is None or BR is None:
        return 0
    mbr_center = (MBR.x + MBR.width / 2, MBR.y + MBR.depth / 2)
    br_center = (BR.x + BR.width / 2, BR.y + BR.depth / 2)
    min_dis = min(MBR.width + BR.width, MBR.depth + BR.depth)
    # diagonal distance is the possible max distance between two points
    max_dis = (boundary.bounds[2] ** 2 + boundary.bounds[3] ** 2) ** 0.5
    distance = (
        (mbr_center[0] - br_center[0]) ** 2 + (mbr_center[1] - br_center[1]) ** 2
    ) ** 0.5

    # normalized_distance = (distance - min_dis) / (max_dis - min_dis)
    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 2. dis(MBR, BA) +
# calculate the distance between the center of the master bedroom and bathroom
# normalize the distance
def dis_MBR_BA(house):
    MBR = None
    BA = None
    for room in house.rooms:
        if room.name == "MBR":
            MBR = room
        if room.name == "BA":
            BA = room
    if MBR is None or BA is None:
        return 0
    mbr_center = (MBR.x + MBR.width / 2, MBR.y + MBR.depth / 2)
    ba_center = (BA.x + BA.width / 2, BA.y + BA.depth / 2)
    min_dis = min(MBR.width + BA.width, MBR.depth + BA.depth)
    max_dis = (boundary.bounds[2] ** 2 + boundary.bounds[3] ** 2) ** 0.5
    distance = (
        (mbr_center[0] - ba_center[0]) ** 2 + (mbr_center[1] - ba_center[1]) ** 2
    ) ** 0.5

    # normalized_distance = (distance - min_dis) / (max_dis - min_dis)
    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 3. orientation of LR +
# (1: north; 0: south) , 1 is preferred (in Southern hemisphere)
# check the orientation of living room
def check_LR_orientation(house):
    LR = None
    for room in house.rooms:
        if room.name == "LR":
            LR = room
    if LR is None:
        return 0

    # check if the living room is in the North (no other rooms are in its north)
    for other_room in house.rooms:
        if other_room.name != "LR" and other_room.name != "BAL":
            if (
                other_room.x < room.x + LR.width
                and other_room.x + other_room.width > LR.x
                and other_room.y + other_room.depth >= LR.y + LR.depth
            ):
                return 0
    return 1


## 4. Light in dining room +
# (0: dark; 1: bright), 1 is preferred
# natural light can come into the dining room through the window
def check_DR_natural_light(house):
    DR = None
    for room in house.rooms:
        if room.name == "DR":
            DR = room
    if DR is None:
        return 0

    # check if there is a window in the dining room
    if len(DR.windows) == 0:
        return 0
    else:
        # check if the window is on the exterior wall
        # this wall is not shared with other rooms and facing outside
        for window in DR.windows:
            for room in house.rooms:
                if room.name != "DR" and room.get_polygon().touches(window.coordinates):
                    intersection_part = room.get_polygon().intersection(
                        window.coordinates
                    )
                    if intersection_part.geom_type in ["LineString", "MultiLineString"]:
                        return 0
    return 1


## 5. Ventilation +
# (Ratio of width to depth of floor plan)
# to be maximized
# Calculate the Ration of width to depth of the whole house
def ventilation(house):
    house_width = house.cluster.bounds[2]
    house_depth = house.cluster.bounds[3]

    min_ratio = 0.2
    max_ratio = 5.0

    ratio = house_width / house_depth
    normalized_ratio = (ratio - min_ratio) / (max_ratio - min_ratio)

    return normalized_ratio


## 6. % of north-facing rooms +
# to be maximized
def north_facing_area(house):
    north_facing_area = 0
    total_area = house.boundary.area
    for room in house.rooms:
        # check if the room is facing south, if no other rooms are in the south
        if room.name != "BAL":
            is_north_facing = True
            for other_room in house.rooms:
                if other_room.name != room.name:
                    if (
                        other_room.y + other_room.depth > room.y + room.depth
                        and other_room.x < room.x + room.width
                        and other_room.x + other_room.width > room.x
                    ):
                        is_north_facing = False
                        break
            if is_north_facing:
                north_facing_area += room.width * room.depth

    per_north_facing = north_facing_area / total_area if total_area != 0 else 0

    return per_north_facing


## 7. % of hall -
# Calculate the % of hallway area in the interior_area
def percentage_hall(house):
    hall_area = 0
    interior_area = 0
    for room in house.rooms:
        if room.name == "HW":
            hall_area += room.width * room.depth
        elif room.name != "BAL":
            interior_area += room.width * room.depth

    per_hall = hall_area / interior_area if interior_area != 0 else 0

    return per_hall


## 8. % of balcony +
# Calculate the % of balcony area in the interior_area
def percentage_balcony(house):
    balcony_area = 0
    interior_area = 0
    for room in house.rooms:
        if room.name == "BAL":
            balcony_area += room.width * room.depth
        elif room.name != "BAL":
            interior_area += room.width * room.depth

    per_balcony = balcony_area / interior_area if interior_area != 0 else 0

    return per_balcony


## 9. Efficiency rate of the house +
# the percentage of interior area to the total area of floor plan
def percentage_interior(house):
    interior_cluster = None
    total_area = house.boundary.area
    for room in house.rooms:
        if room.name != "BAL":
            if interior_cluster is None:
                interior_cluster = room.get_polygon()
            else:
                interior_cluster = interior_cluster.union(room.get_polygon())

    per_interior = interior_cluster.area / total_area

    return per_interior


## 10. DIS (BR, BA) -
# walking distance between the geometric centers of bedroom and bathroom
# to be minimized
def dis_BR_BA(house):
    BR = None
    BA = None
    for room in house.rooms:
        if room.name.startswith("BR"):
            BR = room
        if room.name == "BA":
            BA = room
    if BR is None or BA is None:
        return 0

    # doors of bedroom and bathroom
    if len(BR.doors) > 0 and len(BA.doors) > 0:
        br_door = BR.doors[0]
        ba_door = BA.doors[0]
    else:
        return 0

    br_door_center = br_door.get_center()
    ba_door_center = ba_door.get_center()

    # geometric center of bedroom
    br_center = (BR.x + BR.width / 2, BR.y + BR.depth / 2)
    # geometric center of bathroom
    ba_center = (BA.x + BA.width / 2, BA.y + BA.depth / 2)

    # walking distance: br_center -> br_door_center -> ba_door_center -> ba_center
    walking_distance = (
        (
            (br_center[0] - br_door_center[0]) ** 2
            + (br_center[1] - br_door_center[1]) ** 2
        )
        ** 0.5
        + (
            abs(br_door_center[0] - ba_door_center[0])
            + abs(br_door_center[1] - ba_door_center[1])
        )
        + (
            (ba_door_center[0] - ba_center[0]) ** 2
            + (ba_door_center[1] - ba_center[1]) ** 2
        )
        ** 0.5
    )
    max_dis = boundary.bounds[2] + boundary.bounds[3]
    min_dis = min(BR.width + BA.width, BR.depth + BA.depth)
    normalized_distance = (walking_distance - min_dis) / (max_dis - min_dis)

    return normalized_distance


## 11. DIS (BA, BAL) -
# to be minimised
# walking distance between the geometric centers of bathroom and living room's door
def dis_BA_LR_door(house):
    BA = None
    LR = None
    for room in house.rooms:
        if room.name == "BA":
            BA = room
        if room.name == "LR":
            LR = room
    if BA is None or LR is None:
        return 0

    # doors of bathroom and living room
    if len(BA.doors) > 0 and len(LR.doors) > 0:
        ba_door = BA.doors[0]
        lr_door = LR.doors[0]
    else:
        return 0

    ba_door_center = ba_door.get_center()
    lr_door_center = lr_door.get_center()

    # geometric center of bathroom
    ba_center = (BA.x + BA.width / 2, BA.y + BA.depth / 2)
    # geometric center of living room
    lr_center = (LR.x + LR.width / 2, LR.y + LR.depth / 2)

    # walking distance: ba_center -> ba_door_center -> lr_door_center -> lr_center
    walking_distance = (
        (
            (ba_center[0] - ba_door_center[0]) ** 2
            + (ba_center[1] - ba_door_center[1]) ** 2
        )
        ** 0.5
        + (
            abs(ba_door_center[0] - lr_door_center[0])
            + abs(ba_door_center[1] - lr_door_center[1])
        )
        + (
            (lr_door_center[0] - lr_center[0]) ** 2
            + (lr_door_center[1] - lr_center[1]) ** 2
        )
        ** 0.5
    )
    max_dis = boundary.bounds[2] + boundary.bounds[3]
    min_dis = min(BA.width + LR.width, BA.depth + LR.depth)
    normalized_distance = (walking_distance - min_dis) / (max_dis - min_dis)

    return normalized_distance


## 12. Close(Garage, entry) +
# (0: different side; 1: same side), 1 is preferred
# check if the garage and entry are on the same side
def check_GAR_side(house):
    """
    Check if the garage and entry are on the same side. (0: different side; 1: same side), 1 is preferred
    """
    GAR = None
    for room in house.rooms:
        if room.name == "GAR":
            GAR = room
    if GAR is None:
        return 0

    # check if the garage is on the same side as the entry
    # 1. check which side the entry is on, entry is represented by a line with 2 points, e.g. [(4.5, 0), (5.5, 0)]
    entry_x = entry.xy[0]  # List of x-coordinates: [4.5, 5.5]
    entry_y = entry.xy[1]  # List of y-coordinates: [0, 0]

    # left or right side
    if entry_x[0] == entry_x[1]:
        if GAR.x == entry_x[0] or GAR.x + GAR.width == entry_x[0]:
            return 1
    # top or bottom side
    elif entry_y[0] == entry_y[1]:
        if GAR.y == entry_y[0] or GAR.y + GAR.depth == entry_y[0]:
            return 1
    return 0


## 13. KIC, LR, DR should be in the same cluster
def check_KIC_LR_DR(house):
    """
    KIC, LR, DR should be in the same cluster
    one of them should touch the other two, at least one
    """
    KIC = None
    LR = None
    DR = None
    score = 0
    total = 3
    for room in house.rooms:
        if room.name == "KIC":
            KIC = room
        if room.name == "LR":
            LR = room
        if room.name == "DR":
            DR = room
    if KIC is None and LR is None and DR is None:
        return 0

    if KIC is not None and LR is not None:
        if KIC.get_polygon().touches(LR.get_polygon()):
            score += 1

    if KIC is not None and DR is not None:
        if KIC.get_polygon().touches(DR.get_polygon()):
            score += 1

    if LR is not None and DR is not None:
        if LR.get_polygon().touches(DR.get_polygon()):
            score += 1

    return score / total


## 14. Rooms can touch the entry
def check_no_entry_touch(house):
    """
    Entry is kept for HW | LR | DR. Other rooms cannot touch the entry,
    **(0: unwanted touch, 1: no touch or allowed touch)**
    but can be touched by one point with other rooms
    """
    for room in house.rooms:
        if room.get_polygon().touches(entry) and (
            room.name != "HW" or room.name != "LR" or room.name != "DR"
        ):
            if room.get_polygon().intersection(entry).geom_type in [
                "LineString",
                "MultiLineString",
            ]:
                return 0
    return 1


## 15. Rooms touches the boundary are preferred
def check_boundary_touch(house):
    """
    Rooms touches the boundary are preferred
    """
    score = 0
    for room in house.rooms:
        polygon = room.get_polygon()
        if (
            polygon.bounds[0] == boundary.bounds[0]
            or polygon.bounds[1] == boundary.bounds[1]
            or polygon.bounds[2] == boundary.bounds[2]
            or polygon.bounds[3] == boundary.bounds[3]
        ):
            score += 1
    return score / len(house.rooms)


def cal_overlap_rate(house):
    """
    Calculate the overlap rate of the house
    """
    rooms = house.rooms
    boundary_area = house.boundary.area
    cluster_area = house.cluster.area

    # -1. normalized overlap rate -
    all_rooms_area = sum(room.width * room.depth for room in rooms)
    overlap_rate = (all_rooms_area - cluster_area) / all_rooms_area
    return overlap_rate


def fitness_partial(house):
    """
    Fitness function with overlap penalty. No window and door involved.
    """
    if isinstance(house.cluster, MultiPolygon) or check_GAR_side(house) == 0:
        return 0

    fitness = 0
    # Calculate fitness with special order
    if (fitness := (1 - cal_overlap_rate(house))) < 1:
        return fitness
    elif (fitness := fitness + percentage_interior(house)) < 1.6:
        return fitness
    else:
        fitness += (
            dis_MBR_BR(house) * 2.0
            + dis_MBR_BA(house) * 2.0
            + check_LR_orientation(house)
            # + check_DR_natural_light(house)
            # + ventilation(house)
            + north_facing_area(house)
            # + (1 - percentage_hall(house))
            # + percentage_balcony(house)
            # +percentage_interior(house) * 10.0
            # + (1 - dis_BR_BA(house))
            # + (1 - dis_BA_BAL(house))
            + check_GAR_side(house) * 5.0
            + check_KIC_LR_DR(house)
            + check_no_entry_touch(house)
            # + check_boundary_touch(house) * 5.0
        )

    return fitness


def fitness_all(house):
    """
    Fitness function with overlap penalty. Window and door involved.
    """
    # if isinstance(house.cluster, MultiPolygon) or check_GAR_side(house) == 0:
    if check_GAR_side(house) == 0:
        return 0

    fitness = 0
    # Calculate fitness with special order
    if (fitness := (1 - cal_overlap_rate(house))) < 1:
        return fitness
    elif (fitness := fitness + percentage_interior(house)) < 1.6:
        return fitness
    else:
        fitness += (
            dis_MBR_BR(house) * 2.0
            + dis_MBR_BA(house) * 2.0
            + check_LR_orientation(house)
            + check_DR_natural_light(house)
            + ventilation(house)
            + north_facing_area(house)
            # + (1 - percentage_hall(house))
            # + percentage_balcony(house)
            # + percentage_interior(house) * 10.0
            + (1 - dis_BR_BA(house))
            + (1 - dis_BA_LR_door(house))
            # + check_GAR_side(house) * 5.0
            + check_KIC_LR_DR(house)
            + check_no_entry_touch(house)
            # + check_boundary_touch(house) * 5.0
        )

    return fitness


# Parameters
granularity_size = 0.5
granularity_position = 0.5

particle_num = 200
pso_iter = 200
mcts_iter = 150


# draw windows and doors
def draw_windows_doors(room, ax):
    for window in room.windows:
        x, y = window.coordinates.xy
        ax.plot(x, y, color="blue", linewidth=2, alpha=1.0)

    for door in room.doors:
        x, y = door.coordinates.xy
        ax.plot(x, y, color="brown", linewidth=3, alpha=1.0)


# draw entry, bold and black
def draw_entry(entry, ax):
    x, y = entry.xy
    ax.plot(x, y, color="black", linewidth=5, alpha=1.0)


# Function to save snapshots
def save_snapshots(
    starting_time,
    current_layout,
    iteration_count,
    fitness,
    mcts_iter,
    num_par,
    pso_iter,
):
    plt.figure(figsize=(4, 6))
    x, y = current_layout.boundary.exterior.xy
    plt.plot(x, y, color="black")
    plt.fill(x, y, color="grey", alpha=0.1)
    for i, room in enumerate(current_layout.rooms):
        room_polygon = room.get_polygon()
        x, y = room_polygon.exterior.xy
        plt.fill(x, y, alpha=0.7)
        plt.plot(x, y, alpha=0.9)
        plt.text(
            room_polygon.centroid.x,
            room_polygon.centroid.y,
            room.name,
            ha="center",
            va="center",
        )
        draw_windows_doors(room, plt)
    draw_entry(entry, plt)

    plt.gca().set_aspect("equal", adjustable="box")
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)
    plt.grid(True, linestyle="--", alpha=0.1)
    plt.title(f"Fitness: {fitness}")
    # create target folder if not exists
    if not os.path.exists(f"results/{starting_time}/{mcts_iter}_{num_par}_{pso_iter}"):
        os.makedirs(f"results/{starting_time}/{mcts_iter}_{num_par}_{pso_iter}")
    plt.savefig(
        f"results/{starting_time}/{mcts_iter}_{num_par}_{pso_iter}/{iteration_count}.png"
    )

    plt.close()


# swap the position of random two rooms
def swap_2rooms(layout):
    rooms = layout.rooms
    if len(rooms) < 2:
        return layout  # No need to swap if there are fewer than 2 rooms

    # Filter swappable room pairs that fit within the layout boundaries
    swappable_rooms = []
    for i in range(len(rooms)):
        for j in range(i + 1, len(rooms)):
            room1, room2 = rooms[i], rooms[j]
            if (
                room1.x + room2.width <= layout.boundary.bounds[2]
                and room2.x + room1.width <= layout.boundary.bounds[2]
                and room1.y + room2.depth <= layout.boundary.bounds[3]
                and room2.y + room1.depth <= layout.boundary.bounds[3]
            ):
                swappable_rooms.append((room1, room2))

    if not swappable_rooms:
        return layout  # No swappable room pairs found

    # Randomly select a swappable room pair
    room1, room2 = random.choice(swappable_rooms)

    # Record original positions
    room1_x, room1_y = room1.x, room1.y
    room2_x, room2_y = room2.x, room2.y

    # Swap room positions
    room1.x, room2.x = room2_x, room1_x
    room1.y, room2.y = room2_y, room1_y

    # clear all windows and doors
    for room in rooms:
        room.windows = []
        room.doors = []
    layout.generate_doors_windows()
    layout.cluster = layout.get_cluster()

    # Debugging information
    # print(f" Swapped rooms: {room1.name} and {room2.name}")
    # print(
    #     f" New positions: {room1.name} at ({room1.x}, {room1.y}), {room2.name} at ({room2.x}, {room2.y})"
    # )
    # print(f" Swapped cluster type: {type(layout.cluster)}")

    return layout


# MCTS Node
class MCTSNode:
    def __init__(self, state, parent=None, action=None):
        self.state = state
        self.parent = parent
        self.action = action
        self.children = []
        self.visits = 0
        self.value = 0
        self.untried_actions = state.get_legal_actions()


# Layout State
class LayoutState:
    def __init__(self, boundary, rooms):
        self.boundary = boundary
        self.rooms = rooms
        self.current_room = 0
        self.placed_rooms = []
        self.cluster = None

    def copy(self):
        """Create a deep copy of the state."""
        new_state = LayoutState(self.boundary, self.rooms)
        new_state.current_room = self.current_room
        new_state.placed_rooms = self.placed_rooms.copy()
        new_state.cluster = self.cluster
        return new_state

    def get_legal_actions(self):
        if self.current_room >= len(self.rooms):
            return []

        room_type, width, depth = self.rooms[self.current_room]
        legal_positions = []
        for x in np.arange(0, self.boundary.bounds[2] - width, granularity_position):
            for y in np.arange(
                0, self.boundary.bounds[3] - depth, granularity_position
            ):
                legal_positions.append((x, y))
        # # if cluster is None, all positions are legal
        # if self.cluster is None:
        #     for x in np.arange(0, self.boundary.bounds[2]-width, granularity):
        #         for y in np.arange(0, self.boundary.bounds[3]-depth, granularity):
        #             legal_positions.append((x, y))
        # # if cluster is not None, only positions that not in the cluster are legal
        # else:
        #     for x in np.arange(0, self.boundary.bounds[2]-width, granularity):
        #         for y in np.arange(0, self.boundary.bounds[3]-depth, granularity):
        #             if not Point(x, y).within(self.cluster):
        #                 legal_positions.append((x, y))

        return legal_positions

    def place_room(self, position):
        x, y = position
        room_type, width, depth = self.rooms[self.current_room]
        new_room = Room(room_type, width, depth, x, y)

        # Add buffer to handle geometric issues 💡
        try:
            new_poly = new_room.get_polygon().buffer(0)
            if self.cluster:
                combined = self.cluster.buffer(0).union(new_poly)
                if not combined.is_valid:
                    return False
                self.cluster = combined
            else:
                self.cluster = new_poly

            self.placed_rooms.append(new_room)
            # print("Placed room:", len(self.placed_rooms))
            self.current_room += 1
            # print("Current room:", self.current_room)
            return True
        except Exception:
            return False

    def all_placed(self):
        """Check if all rooms have been placed."""
        return self.current_room >= len(self.rooms)

    def get_house(self):
        """Create a House object from the current state."""
        # print("Placed--rooms:", len(self.placed_rooms))
        return House(self.placed_rooms, self.boundary)


# MCTS algorithm
class MCTS:
    def __init__(self, boundary, rooms, iterations=mcts_iter, exploration=1.414):
        self.boundary = boundary
        self.rooms = rooms
        self.iterations = iterations
        self.exploration = exploration

    def search(self):
        root_state = LayoutState(self.boundary, self.rooms)
        root_node = MCTSNode(root_state)
        best_reward = -float("inf")
        best_state = None

        for _ in range(self.iterations):
            node = root_node
            state = root_state.copy()
            valid_simulation = True

            # Selection with validity check
            while node.untried_actions == [] and node.children != []:
                # print("Selecting child")
                node = self.select_child(node)
                if not state.place_room(node.action):
                    valid_simulation = False
                    break

            if not valid_simulation:
                continue  # Skip invalid paths

            # Expansion with smart action selection
            if node.untried_actions:
                action = random.choice(node.untried_actions)
                # print("Expanding with action")
                if state.place_room(action):
                    node.untried_actions.remove(action)
                    new_node = MCTSNode(state.copy(), node, action)
                    node.children.append(new_node)
                    node = new_node

            # Enhanced simulation with retry mechanism
            simulation_state = state.copy()
            max_attempts = 100
            attempts = 0
            while not simulation_state.all_placed() and attempts < max_attempts:
                legal_actions = simulation_state.get_legal_actions()
                if not legal_actions:
                    valid_simulation = False
                    break
                action = random.choice(legal_actions)
                if not simulation_state.place_room(action):
                    valid_simulation = False
                    break
                # print("Correct place?:", simulation_state.place_room(action))

                attempts += 1

            # Reward calculation
            if not valid_simulation:
                reward = 0.0
            elif simulation_state.all_placed():
                house = simulation_state.get_house()
                house.generate_doors_windows()  # Generate windows and doors for fitness calculation
                reward = fitness_all(house)

                # Persist the generated windows and doors in the simulation state
                simulation_state.placed_rooms = (
                    house.rooms
                )  # Update rooms with windows and doors
                simulation_state.cluster = house.cluster  # Update the cluster

                if reward > best_reward:
                    best_reward = reward
                    best_state = simulation_state

            # Backpropagation
            current_node = node
            while current_node is not None:
                current_node.visits += 1
                current_node.value += reward
                current_node = current_node.parent

        # print("Reward:", reward)
        # # print("node.visits:", root_node.visits)
        # print("node.children:", len(root_node.children))
        # # print("node.action:", node.action)
        # return self.get_best_layout(root_node.children)
        if best_state:
            return best_state.get_house()
        else:
            return self.get_best_layout(root_node.children)

    def select_child(self, node):
        log_total = math.log(node.visits)
        best_score = -float("inf")
        best_child = None

        for child in node.children:
            score = child.value / child.visits + self.exploration * math.sqrt(
                log_total / child.visits
            )
            if score > best_score:
                best_score = score
                best_child = child

        return best_child

    def get_best_layout(self, children):
        best_value = -float("inf")
        best_node = None
        for child in children:
            if child.visits > best_value:
                best_value = child.visits
                best_node = child
        return best_node.state.get_house()


# PSO-MCTS algorithm
class PSO:
    def __init__(self, rooms_range, num_particles, max_iter, mcts_iter):
        self.rooms_range = rooms_range
        self.num_particles = num_particles
        self.max_iter = max_iter
        self.gbest_fitness = -float("inf")
        self.gbest_sizes = []
        self.gbest_layout = None
        self.mcts_iter = mcts_iter

        # Initialize particles
        self.particles = []
        for _ in range(num_particles):
            particle = {
                "sizes": [self.random_size(room) for room in rooms_range],
                "velocity": [[0, 0] for _ in rooms_range],
                "pbest_fitness": -float("inf"),
                # "pbest_fitness_partial": -float("inf"),
                "pbest_sizes": [],
            }
            self.particles.append(particle)

    def random_size(self, room):
        min_w, max_w, min_d, max_d = room[1:]
        return [
            round(random.uniform(min_w, max_w) / granularity_size) * granularity_size,
            round(random.uniform(min_d, max_d) / granularity_size) * granularity_size,
        ]

    def optimize(self, boundary):
        current_time = datetime.now().strftime("%Y%m%d-%H%M%S")
        for _ in range(self.max_iter):
            for particle in self.particles:
                # Evaluate fitness using MCTS
                rooms = []
                for i, size in enumerate(particle["sizes"]):
                    rooms.append((self.rooms_range[i][0], size[0], size[1]))

                mcts = MCTS(boundary, rooms)
                best_layout = mcts.search()
                # print("Best layout:", len(best_layout.rooms))
                fitness = fitness_all(best_layout)

                # Update personal best
                if fitness > particle["pbest_fitness"]:
                    particle["pbest_fitness"] = fitness
                    particle["pbest_sizes"] = particle["sizes"].copy()

                # Update global best
                if fitness > self.gbest_fitness:
                    self.gbest_fitness = fitness
                    self.gbest_sizes = particle["sizes"].copy()
                    self.gbest_layout = best_layout

            # local search: swap the position of random two rooms
            swaped_layout = swap_2rooms(copy.deepcopy(self.gbest_layout))
            swaped_fitness = fitness_all(swaped_layout)
            if swaped_fitness > self.gbest_fitness:
                self.gbest_layout = swaped_layout
                self.gbest_fitness = swaped_fitness

            # Update velocities and positions
            for particle in self.particles:
                for i in range(len(particle["sizes"])):
                    w = 0.5  # inertia
                    c1 = 0.6  # cognitive
                    c2 = 0.9  # social

                    v = [
                        w * particle["velocity"][i][0]
                        + c1
                        * random.random()
                        * (particle["pbest_sizes"][i][0] - particle["sizes"][i][0])
                        + c2
                        * random.random()
                        * (self.gbest_sizes[i][0] - particle["sizes"][i][0]),
                        w * particle["velocity"][i][1]
                        + c1
                        * random.random()
                        * (particle["pbest_sizes"][i][1] - particle["sizes"][i][1])
                        + c2
                        * random.random()
                        * (self.gbest_sizes[i][1] - particle["sizes"][i][1]),
                    ]

                    particle["velocity"][i] = v
                    particle["sizes"][i] = [
                        max(
                            min(
                                round(
                                    (particle["sizes"][i][0] + v[0]) / granularity_size
                                )
                                * granularity_size,
                                self.rooms_range[i][2],
                            ),
                            self.rooms_range[i][1],
                        ),
                        max(
                            min(
                                round(
                                    (particle["sizes"][i][1] + v[1]) / granularity_size
                                )
                                * granularity_size,
                                self.rooms_range[i][4],
                            ),
                            self.rooms_range[i][3],
                        ),
                    ]

            # Save snapshots, the initial, final and every 25 iterations
            if _ == 0 or _ == self.max_iter - 1 or _ % 25 == 0:
                # force the count number to be like: 0000, 0025, 0050, ...
                save_snapshots(
                    current_time,
                    self.gbest_layout,
                    str(_).zfill(4),
                    self.gbest_fitness,
                    self.mcts_iter,
                    self.num_particles,
                    self.max_iter,
                )
                print(
                    f"Iteration: {_}, Polygon Type: {self.gbest_layout.cluster.geom_type}, Fitness: {self.gbest_fitness}"
                )

        return self.gbest_layout


# List of configurations
# (mcts_iter, num_particles, pso_iter)
configurations = [
    # (1000, 200, 200),
    # (150, 200, 200),
    # (150, 200, 2000),
    # (150, 600, 200),
    # (300, 20, 2000),
    # (300, 200, 200),
    # (300, 200, 2000),
    # (300, 200, 200),
    (600, 400, 400),
]


def run_pso_mcts(mcts_iter, num_particles, pso_iter):
    print(
        f"Running PSO-MCTS with mcts_iter={mcts_iter}, num_particles={num_particles}, pso_iter={pso_iter}"
    )
    # algorithm logic
    pso = PSO(
        rooms_range, num_particles=num_particles, max_iter=pso_iter, mcts_iter=mcts_iter
    )
    best_house = pso.optimize(boundary)

    # final check
    if len(best_house.rooms) < len(rooms_range):
        print("Failed to place all rooms")

    return best_house


# Function to run a single configuration
def run_config(config):
    mcts_iter, particles, pso_iter = config
    run_pso_mcts(mcts_iter, particles, pso_iter)


# Use multiprocessing to run configurations in parallel
if __name__ == "__main__":
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(run_config, configurations)
