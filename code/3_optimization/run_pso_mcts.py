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
house_width = 16
house_depth = 10

boundary = Polygon(
    [(0, 0), (house_width, 0), (house_width, house_depth), (0, house_depth)]
)

# The entry position
entry = LineString([(0, 5), (0, 6)])


# Define room types (type, min_width, max_width, min_depth, max_depth)
rooms_range = [
    ("GAR", 5, 6, 3, 6),  # Garage
    ("LDR", 2, 3, 2, 3),  # Laundry Room
    ("DR", 3, 5, 3, 5),  # Dining Room
    ("LR", 3, 6, 3, 6),  # Living Room
    ("KIC", 2, 4, 2, 4),  # Kitchen
    ("MBR", 3, 6, 3, 6),  # Master Bedroom
    ("BR1", 3, 5, 3, 5),  # Bedroom
    ("BA", 2, 4, 2, 4),  # Bathroom
    # ("BR2", 3, 5, 3, 5),  # Bedroom2
    # ("HW", 2, 10, 1, 2),  # Hallway
    # ("BAL", 2, 8, 1, 2),  # Balcony
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
        self.width = int(width)  # Ensure width is an integer
        self.depth = int(depth)  # Ensure depth is an integer
        self.x = int(x)  # Ensure x-coordinate is an integer
        self.y = int(y)  # Ensure y-coordinate is an integer
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


# ====================Fitness functions=========================================


## 1. dis(MBR, BR) +
def dis_MBR_BR(house):
    """
    1. dis(MBR, BR) +\\
    calculate the distance between the center of the master bedroom and other bedrooms. Normalized
    """
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
    min_dis = min((MBR.width + BR.width) / 2, (MBR.depth + BR.depth) / 2)
    # 2 rooms put at 2 opposite corners of the house
    max_dis = (
        (boundary.bounds[2] - (MBR.width + BR.width) / 2) ** 2
        + (boundary.bounds[3] - (MBR.depth + BR.depth) / 2) ** 2
    ) ** 0.5
    distance = (
        (mbr_center[0] - br_center[0]) ** 2 + (mbr_center[1] - br_center[1]) ** 2
    ) ** 0.5

    normalized_distance = (distance - min_dis) / (max_dis - min_dis)

    return normalized_distance


## 2. dis(MBR, BA) +
def dis_MBR_BA(house):
    """
    2. dis(MBR, BA) + \\
    calculate the distance between the center of the master bedroom and bathroom. Normalized
    """
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
    min_dis = min((MBR.width + BA.width) / 2, (MBR.depth + BA.depth) / 2)
    # 2 rooms put at 2 opposite corners of the house
    max_dis = (
        (boundary.bounds[2] - (MBR.width + BA.width) / 2) ** 2
        + (boundary.bounds[3] - (MBR.depth + BA.depth) / 2) ** 2
    ) ** 0.5
    distance = (
        (mbr_center[0] - ba_center[0]) ** 2 + (mbr_center[1] - ba_center[1]) ** 2
    ) ** 0.5

    normalized_distance = (distance - min_dis) / (max_dis - min_dis)

    return normalized_distance


## 3. orientation of LR +
def check_LR_orientation(house):
    """
    3. orientation of LR +\\
    check the orientation of living room.
    (1: north; 0: south) , 1 is preferred (in Southern hemisphere). Normalized.
    """
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
def check_DR_natural_light(house):
    """
    4. Light in dining room +\\
    check if there is a window in the dining room so that natural light can come in.
    (0: dark; 1: bright), 1 is preferred
    """
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
def ventilation(house):
    """
    5. Ventilation +\\
    (Ratio of width to depth of floor plan)\\
    to be maximized\\
    Calculate the Ration of width to depth of the whole house
    """
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
    """
    6. % of north-facing rooms +\\
    Calculate the % of north-facing rooms area in the interior_area
    """
    north_facing_area = 0
    total_area = 0
    for room in house.rooms:
        # check if the room is facing south, if no other rooms are in the south
        if room.name != "BAL":
            room_area = room.width * room.depth
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
            total_area += room_area

            if is_north_facing:
                north_facing_area += room_area

    per_north_facing = north_facing_area / total_area if total_area != 0 else 0

    return per_north_facing


## 7. % of hall -
def percentage_hall(house):
    """
    7. % of hall -\\
    Calculate the % of hallway area in the interior_area
    """
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
def percentage_balcony(house):
    """
    8. % of balcony +\\
    Calculate the % of balcony area in the interior_area
    """
    balcony_area = 0
    interior_area = 0
    for room in house.rooms:
        if room.name == "BAL":
            balcony_area += room.width * room.depth
        elif room.name != "BAL":
            interior_area += room.width * room.depth

    per_balcony = balcony_area / interior_area if interior_area != 0 else 0

    return per_balcony


## 9. Utilization rate of the house +
def utilization_rate(house):
    """
    9. Utilization rate of the house +\\
    Calculate the % of interior area in the total area of the floor plan
    """
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
def dis_BR_BA(house):
    """
    10. DIS (BR, BA) -\\
    walking distance between the geometric centers of bedroom and bathroom
    """
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


## 11. DIS (BA, LR_door) -
def dis_BA_LR_door(house):
    """
    11. DIS (BA, LR_door) -\\
    walking distance between the geometric centers of bathroom and living room's door
    """
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


## 13.1 dis(KIC, LR) -
def dis_KIC_LR(house):
    """
    13.1 dis(KIC, LR) -\\
    Calculate the normalized distance between the Kitchen (KIC) and Living Room (LR). Closer is better
    """
    KIC = None
    LR = None
    for room in house.rooms:
        if room.name == "KIC":
            KIC = room
        if room.name == "LR":
            LR = room
    if KIC is None or LR is None:
        return 0

    KIC_polygon = KIC.get_polygon()
    LR_polygon = LR.get_polygon()
    if (
        KIC_polygon.within(LR_polygon)
        or LR_polygon.within(KIC_polygon)
        or KIC_polygon.overlaps(LR_polygon)
    ):
        return 1

    distance = KIC_polygon.distance(LR_polygon)
    max_dis = (
        (boundary.bounds[2] - (KIC.width + LR.width)) ** 2
        + (boundary.bounds[3] - (KIC.depth + LR.depth)) ** 2
    ) ** 0.5

    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 13.2 dis(KIC, DR) -
def dis_KIC_DR(house):
    """
    13.2 dis(KIC, DR) -\\
    Calculate the normalized distance between the Kitchen (KIC) and Dining Room (DR). Closer is better
    """
    KIC = None
    DR = None
    for room in house.rooms:
        if room.name == "KIC":
            KIC = room
        if room.name == "DR":
            DR = room
    if KIC is None or DR is None:
        return 0

    KIC_polygon = KIC.get_polygon()
    DR_polygon = DR.get_polygon()
    if (
        KIC_polygon.within(DR_polygon)
        or DR_polygon.within(KIC_polygon)
        or KIC_polygon.overlaps(DR_polygon)
    ):
        return 1

    distance = KIC_polygon.distance(DR_polygon)
    max_dis = (
        (boundary.bounds[2] - (KIC.width + DR.width)) ** 2
        + (boundary.bounds[3] - (KIC.depth + DR.depth)) ** 2
    ) ** 0.5

    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 13.3 dis(LR, DR) -
def dis_LR_DR(house):
    """
    13.3 dis(LR, DR) -\\
    Calculate the normalized distance between the Living Room (LR) and Dining Room (DR). Closer is better
    """
    LR = None
    DR = None
    for room in house.rooms:
        if room.name == "LR":
            LR = room
        if room.name == "DR":
            DR = room
    if LR is None or DR is None:
        return 0

    LR_polygon = LR.get_polygon()
    DR_polygon = DR.get_polygon()
    if (
        LR_polygon.within(DR_polygon)
        or DR_polygon.within(LR_polygon)
        or LR_polygon.overlaps(DR_polygon)
    ):
        return 1

    distance = LR_polygon.distance(DR_polygon)
    max_dis = (
        (boundary.bounds[2] - (LR.width + DR.width)) ** 2
        + (boundary.bounds[3] - (LR.depth + DR.depth)) ** 2
    ) ** 0.5

    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 14. Rooms can touch the entry
def check_no_entry_touch(house):
    """
    Entry is kept for HW, LR, DR, KIC. Other rooms cannot touch the entry,
    **(0: unwanted touch, 1: no touch or allowed touch)**
    but can be touched by one point with other rooms
    """
    for room in house.rooms:
        if room.get_polygon().touches(entry) and (
            room.name != "HW"
            and room.name != "LR"
            and room.name != "DR"
            and room.name != "KIC"
        ):
            if room.get_polygon().intersection(entry).geom_type in [
                "LineString",
                "MultiLineString",
            ]:
                return 0
    return 1


## 15. Calculate the overlap rate of the layout
def cal_overlap_rate(house):
    """
    15. cal_overlap_rate - \\
    Calculate the overlap rate of the house. To be minimized
    """
    rooms = house.rooms
    boundary_area = house.boundary.area
    cluster_area = house.cluster.area

    # 16. normalized overlap rate -
    all_rooms_area = sum(room.width * room.depth for room in rooms)
    overlap_rate = (all_rooms_area - cluster_area) / all_rooms_area
    return overlap_rate


## 16. Differentiate public and private rooms
def diff_public_private(house):
    """
    16. Differentiate public and private rooms
    """
    public_rooms = ["LR", "DR", "KIC", "LDR", "GAR"]
    private_rooms = ["MBR", "BA", "BR1"]
    total = len(public_rooms) * len(private_rooms)
    score = 0

    # The center of the entry is the middle of the line
    entry_center = entry.interpolate(0.5, normalized=True)
    for room in house.rooms:
        if room.name in public_rooms:
            for other_room in house.rooms:
                if other_room.name in private_rooms:
                    # calculate the distance between the entry and 2 rooms
                    room_center = (room.x + room.width / 2, room.y + room.depth / 2)
                    other_room_center = (
                        other_room.x + other_room.width / 2,
                        other_room.y + other_room.depth / 2,
                    )
                    # distance between entry and public room
                    dis_entry_room = (entry_center.x - room_center[0]) ** 2 + (
                        entry_center.y - room_center[1]
                    ) ** 2
                    # distance between entry and private room
                    dis_entry_other_room = (
                        entry_center.x - other_room_center[0]
                    ) ** 2 + (entry_center.y - other_room_center[1]) ** 2
                    # compare the distance
                    if dis_entry_room < dis_entry_other_room:
                        score += 1

    # Normalize the score
    score = score / total

    return score


## 17. KIC facing north +
def check_KIC_orientation(house):
    """
    17. KIC facing north +\\
    check if the kitchen is in the North (no other rooms are in its north).
    (1: north; 0: south) , 1 is preferred (in Southern hemisphere). Normalized.
    """
    KIC = None
    for room in house.rooms:
        if room.name == "KIC":
            KIC = room
    if KIC is None:
        return 0
    # check if the kitchen is in the North (no other rooms are in its north)
    for other_room in house.rooms:
        if other_room.name != "KIC" and other_room.name != "BAL":
            if (
                other_room.x < room.x + KIC.width
                and other_room.x + other_room.width > KIC.x
                and other_room.y + other_room.depth >= KIC.y + KIC.depth
            ):
                return 0
    return 1


## 18. dis(LDR, KIC) -
def dis_LDR_KIC(house):
    """
    18. dis(LDR, KIC) -\\
    calculate the distance between the center of the LDR and KIC. Normalized
    """
    LDR = None
    KIC = None
    for room in house.rooms:
        if room.name == "LDR":
            LDR = room
        if room.name == "KIC":
            KIC = room
    if LDR is None or KIC is None:
        return 0

    LDR_polygon = LDR.get_polygon()
    KIC_polygon = KIC.get_polygon()
    if (
        LDR_polygon.within(KIC_polygon)
        or KIC_polygon.within(LDR_polygon)
        or LDR_polygon.overlaps(KIC_polygon)
    ):
        return 1

    distance = LDR_polygon.distance(KIC_polygon)
    # 2 rooms put at 2 opposite corners of the house
    max_dis = (
        (boundary.bounds[2] - (LDR.width + KIC.width)) ** 2
        + (boundary.bounds[3] - (LDR.depth + KIC.depth)) ** 2
    ) ** 0.5

    # normalized_distance = (distance - min_dis) / (max_dis - min_dis)
    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 19. dis(Garage, KIC) -
def dis_GAR_KIC(house):
    """
    19. dis(GAR, KIC) -\\
    calculate the distance between the GAR and KIC. Normalized
    """
    GAR = None
    KIC = None
    for room in house.rooms:
        if room.name == "GAR":
            GAR = room
        if room.name == "KIC":
            KIC = room
    if GAR is None or KIC is None:
        return 0

    GAR_polygon = GAR.get_polygon()
    KIC_polygon = KIC.get_polygon()
    if (
        GAR_polygon.within(KIC_polygon)
        or KIC_polygon.within(GAR_polygon)
        or GAR_polygon.overlaps(KIC_polygon)
    ):
        return 1

    distance = GAR_polygon.distance(KIC_polygon)
    # 2 rooms put at 2 opposite corners of the house
    max_dis = (
        (boundary.bounds[2] - (GAR.width + KIC.width)) ** 2
        + (boundary.bounds[3] - (GAR.depth + KIC.depth)) ** 2
    ) ** 0.5

    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 20. dis(MBR, LR) +
def dis_MBR_LR(house):
    """
    20. dis(MBR, LR) +\\
    calculate the distance between the center of the master bedroom and living room. Normalized
    """
    MBR = None
    LR = None
    for room in house.rooms:
        if room.name == "MBR":
            MBR = room
        if room.name == "LR":
            LR = room
    if MBR is None or LR is None:
        return 0

    MBR_polygon = MBR.get_polygon()
    LR_polygon = LR.get_polygon()
    if (
        MBR_polygon.within(LR_polygon)
        or LR_polygon.within(MBR_polygon)
        or MBR_polygon.overlaps(LR_polygon)
    ):
        return 0

    distance = MBR_polygon.distance(LR_polygon)
    # 2 rooms put at 2 opposite corners of the house
    max_dis = (
        (boundary.bounds[2] - (MBR.width + LR.width)) ** 2
        + (boundary.bounds[3] - (MBR.depth + LR.depth)) ** 2
    ) ** 0.5
    # normalized_distance = (distance - min_dis) / (max_dis - min_dis)
    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


## 21. dis(MBR, KIC) +
def dis_MBR_KIC(house):
    """
    21. dis(MBR, KIC) +\\
    calculate the distance between the center of the master bedroom and kitchen. Normalized
    """
    MBR = None
    KIC = None
    for room in house.rooms:
        if room.name == "MBR":
            MBR = room
        if room.name == "KIC":
            KIC = room
    if MBR is None or KIC is None:
        return 0

    MBR_polygon = MBR.get_polygon()
    KIC_polygon = KIC.get_polygon()
    if (
        MBR_polygon.within(KIC_polygon)
        or KIC_polygon.within(MBR_polygon)
        or MBR_polygon.overlaps(KIC_polygon)
    ):
        return 0

    distance = MBR_polygon.distance(KIC_polygon)
    # 2 rooms put at 2 opposite corners of the house
    max_dis = (
        (boundary.bounds[2] - (MBR.width + KIC.width)) ** 2
        + (boundary.bounds[3] - (MBR.depth + KIC.depth)) ** 2
    ) ** 0.5
    # normalized_distance = (distance - min_dis) / (max_dis - min_dis)
    normalized_distance = (distance - 0) / (max_dis - 0)

    return normalized_distance


def fitness_partial(house):
    """
    Partial Fitness function with overlap penalty & utilization rate. No window and door involved.
    """
    fitness = 0
    penalty = 0
    if check_GAR_side(house) == 0:
        penalty += 0.5
    if check_no_entry_touch(house) == 0:
        penalty += 0.5

    # Calculate fitness with special order
    if (fitness := (1 - cal_overlap_rate(house))) < 1:
        return fitness - penalty
    elif (fitness := fitness + utilization_rate(house)) < 1.75:
        return fitness - penalty
    else:
        fitness += (
            # utilization_rate(house) * 20
            # + dis_MBR_BR(house)
            # + dis_MBR_BA(house)
            check_LR_orientation(house) * 2
            # # + ventilation(house)
            # + north_facing_area(house)
            # # + check_GAR_side(house)
            + (1 - dis_KIC_LR(house))
            + (1 - dis_KIC_DR(house))
            + (1 - dis_LR_DR(house))
            # # + check_no_entry_touch(house)
            + diff_public_private(house) * 15
            + check_KIC_orientation(house) * 2
            + (1 - dis_LDR_KIC(house))
            + (1 - dis_GAR_KIC(house))
            + dis_MBR_LR(house)
            + dis_MBR_KIC(house)
        )

    # fitness += diff_public_private(house)

    return fitness - penalty


def fitness_all(house):
    """
    Fitness function with overlap penalty. Windows and doors involved.
    """
    # if isinstance(house.cluster, MultiPolygon) or check_GAR_side(house) == 0:
    if check_GAR_side(house) == 0 or check_no_entry_touch(house) == 0:
        return 0

    fitness = 0
    # Calculate fitness with special order
    if (fitness := (1 - cal_overlap_rate(house))) < 1:
        return fitness
    # elif (fitness := fitness + utilization_rate(house)) < 1.8:
    #     return fitness
    else:
        fitness += (
            utilization_rate(house)
            + dis_MBR_BR(house)
            # + dis_MBR_BA(house)
            # + check_LR_orientation(house)
            # # + ventilation(house)
            # + north_facing_area(house)
            # # + check_GAR_side(house)
            # + (1 - dis_KIC_LR(house))
            # + (1 - dis_KIC_DR(house))
            # + (1 - dis_LR_DR(house))
            # # + check_no_entry_touch(house)
            # + diff_public_private(house) * 12
            # + check_KIC_orientation(house)
            # + (1 - dis_LDR_KIC(house))
            # + (1 - dis_GAR_KIC(house))
            # + dis_MBR_LR(house)
            # + dis_MBR_KIC(house)
            # ============⬇️window/door involved⬇️===========================
            ## + check_DR_natural_light(house)
            # + (1 - dis_BR_BA(house))
            # + (1 - dis_BA_LR_door(house))
        )

    return fitness


# Parameters
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
    plt.figure(figsize=(6, 4))
    x, y = current_layout.boundary.exterior.xy
    plt.plot(x, y, color="black")
    plt.fill(x, y, color="grey", alpha=0.1)
    for i, room in enumerate(current_layout.rooms):
        room_polygon = room.get_polygon()
        x, y = room_polygon.exterior.xy
        plt.fill(x, y, alpha=0.7)
        # use dotted line for LR, KIC, DR
        if room.name in ["LR", "KIC", "DR"]:
            plt.plot(x, y, color="black", linestyle="--", alpha=0.8)
        else:
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
    plt.grid(True, linestyle="--", alpha=0.2)
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

    # # clear all windows and doors
    # for room in rooms:
    #     room.windows = []
    #     room.doors = []
    # layout.generate_doors_windows()
    # layout.cluster = layout.get_cluster()

    # Debugging information
    # print(f" Swapped rooms: {room1.name} and {room2.name}")
    # print(
    #     f" New positions: {room1.name} at ({room1.x}, {room1.y}), {room2.name} at ({room2.x}, {room2.y})"
    # )
    # print(f" Swapped cluster type: {type(layout.cluster)}")

    return layout


# 1+1 EA algorithm
class Position_OnePlusOneEA:
    """
    1+1 EA algorithm for mutating room's position
    """

    def __init__(self, boundary, rooms, max_iter=100):
        self.boundary = boundary
        self.rooms = rooms
        self.max_iter = max_iter
        self.best_layout = None
        self.best_fitness = -float("inf")

    def search(self):
        # Initialize the first layout
        placed_rooms = []
        legal_positions_array = []
        index_position_array = []
        for i, room in enumerate(self.rooms):
            # print("Room:", room)

            room_type, width, depth = room
            # Get all possible positions for rooms
            x_positions = range(0, int(self.boundary.bounds[2] - width) + 1)
            y_positions = range(0, int(self.boundary.bounds[3] - depth) + 1)
            legal_positions = [(x, y) for x in x_positions for y in y_positions]
            legal_positions_array.append(legal_positions)
            index_position_array.append(i)
            len_legal_positions = len(legal_positions)
            if len_legal_positions == 0:
                print("No legal positions available for room:", room_type)
                continue

            # Randomly select a position for the room
            index = random.randint(0, len_legal_positions - 1)
            x, y = legal_positions[index]
            room = Room(room_type, width, depth, x, y)
            placed_rooms.append(room)

        # Create a new layout
        new_layout = House(placed_rooms, self.boundary)
        self.best_layout = new_layout
        self.best_fitness = fitness_partial(new_layout)
        # print("Initial fitness:", self.best_fitness)
        # Iterate to find the best layout
        for j in range(self.max_iter):
            for i, room in enumerate(placed_rooms):
                # Generate a mutation rate [0, 1]
                mutation_rate = random.uniform(0, 1)
                # Mutate this room if mutation_rate < 1/8
                if mutation_rate <= 0.125:
                    # Get the room type, width, and depth
                    room_type = room.name
                    width = room.width
                    depth = room.depth
                    # Get all possible positions for the room
                    legal_positions = legal_positions_array[i]
                    len_legal_positions = len(legal_positions)
                    if len_legal_positions == 0:
                        print("No legal positions available for room:", room_type)
                        continue
                    # Get current index of position array
                    index = index_position_array[i]
                    # Generate a random number [-50, 50]
                    random_number = random.randint(-50, 50)
                    # Get the new index of position array
                    new_index = (index + random_number) % len_legal_positions
                    # Get the new position of the room
                    new_position = legal_positions[new_index]
                    # Create a new room with the new position
                    new_room = Room(room_type, width, depth, *new_position)
                    # Generate a new layout with the new room
                    new_layout = House(
                        placed_rooms[:i] + [new_room] + placed_rooms[i + 1 :],
                        self.boundary,
                    )

                    # Calculate the fitness of the new layout
                    new_fitness = fitness_partial(new_layout)
                    # If the new fitness is better than the best fitness, update the best layout
                    if new_fitness > self.best_fitness:
                        self.best_layout = new_layout
                        self.best_fitness = new_fitness
                        placed_rooms = (
                            placed_rooms[:i] + [new_room] + placed_rooms[i + 1 :]
                        )
                        index_position_array[i] = new_index

        return self.best_layout, self.best_fitness


# PSO & (1+1)_EA algorithm
class PSO:
    def __init__(self, rooms_range, num_particles, max_iter, mcts_iter):
        self.rooms_range = rooms_range
        self.num_particles = num_particles
        self.max_iter = max_iter
        self.gbest_fitness = -float("inf")
        self.gbest_sizes = []
        self.gbest_layout = None
        self.mcts_iter = mcts_iter
        self.full_fitness = -float("inf")

        # Initialize particles
        self.particles = []
        for _ in range(num_particles):
            particle = {
                "sizes": [self.random_size(room) for room in rooms_range],
                "velocity": [[0, 0] for _ in rooms_range],
                "pbest_fitness": -float("inf"),
                "pbest_sizes": [],
            }
            self.particles.append(particle)

    def random_size(self, room):
        min_w, max_w, min_d, max_d = room[1:]
        return [
            random.randint(int(min_w), int(max_w)),  # Ensure width is an integer
            random.randint(int(min_d), int(max_d)),  # Ensure depth is an integer
        ]

    def optimize(self, boundary):
        current_time = datetime.now().strftime("%Y%m%d-%H%M%S")
        print(f"PSO started at {current_time}")
        for _ in range(self.max_iter):
            # Evaluate fitness for each particle
            for index, particle in enumerate(self.particles):
                # Evaluate fitness using MCTS
                rooms = []
                for i, size in enumerate(particle["sizes"]):
                    rooms.append((self.rooms_range[i][0], size[0], size[1]))

                # mcts = MCTS(boundary, rooms)
                # best_layout = mcts.search()
                one_one_EA = Position_OnePlusOneEA(boundary, rooms)
                best_layout, fitness = one_one_EA.search()
                # if _ % 5 == 0:
                #     print(f"Particle {index}, Iteration {_}: Best fitness: {fitness}")

                # fitness = fitness_all(best_layout)

                # Update personal best
                if fitness > particle["pbest_fitness"]:
                    particle["pbest_fitness"] = fitness
                    particle["pbest_sizes"] = particle["sizes"].copy()

                # Update global best
                if fitness > self.gbest_fitness:
                    self.gbest_fitness = fitness
                    self.gbest_sizes = particle["sizes"].copy()
                    self.gbest_layout = best_layout

            # Local search: swap the position of random two rooms
            swaped_layout = swap_2rooms(copy.deepcopy(self.gbest_layout))
            swaped_fitness = fitness_partial(swaped_layout)
            if swaped_fitness > self.gbest_fitness:
                self.gbest_layout = swaped_layout
                self.gbest_fitness = swaped_fitness

            # if self.gbest_fitness >= 1.6:
            #     # Add doors and windows to the global best layout
            #     copied_layout = copy.deepcopy(self.gbest_layout)
            #     copied_layout.generate_doors_windows()
            #     copied_layout.get_cluster()
            #     tmp_full_fitness = fitness_all(copied_layout)
            #     if tmp_full_fitness > self.full_fitness:
            #         self.full_fitness = tmp_full_fitness
            #         self.gbest_layout = copied_layout

            print(
                f"Iteration {_} finished at: {datetime.now().strftime('%Y%m%d-%H%M%S')}, Fitness: {self.gbest_fitness}"
            )
            # Update velocities and positions
            for particle in self.particles:
                for i in range(len(particle["sizes"])):
                    w = 0.5  # inertia
                    c1 = 0.4  # cognitive
                    c2 = 0.4  # social

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
                                particle["sizes"][i][0] + int(v[0]),
                                self.rooms_range[i][2],
                            ),
                            self.rooms_range[i][1],
                        ),
                        max(
                            min(
                                particle["sizes"][i][1] + int(v[1]),
                                self.rooms_range[i][4],
                            ),
                            self.rooms_range[i][3],
                        ),
                    ]

            # Save snapshots, the initial, final, and every 25 iterations
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


# (1+1)_EA algorithm for mutate room's width and depth
class Size_OnePlusOneEA:
    """
    1+1 EA algorithm for mutating room's width and depth after 1+1 EA for Position
    """

    def __init__(self, boundary, rooms_range, max_iter=20000):
        self.boundary = boundary
        self.rooms_range = rooms_range
        self.max_iter = max_iter
        self.best_layout = None
        self.best_layout_rooms = []
        self.best_fitness = -float("inf")

    def search(self):
        current_time = datetime.now().strftime("%Y%m%d-%H%M%S")
        # 0. generate the first layout
        rooms = []
        for i, room in enumerate(self.rooms_range):
            w = random.randint(int(room[1]), int(room[2]))
            d = random.randint(int(room[3]), int(room[4]))
            self.best_layout_rooms.append((room[0], w, d))

        # 1. randomly mutate room(s)' width and depth
        for i in range(self.max_iter):
            rooms = self.best_layout_rooms.copy()
            for j, room in enumerate(rooms):
                # generate a mutation rate [0, 1]
                mutation_rate = random.uniform(0, 1)
                # Mutate this room size if mutation_rate < 1/8
                if mutation_rate <= 0.125:
                    # Get the room type, width, and depth
                    room_type = room[0]
                    width = room[1]
                    depth = room[2]
                    # generate 2 random numbers [-3, 3]
                    random_number1 = random.randint(-3, 3)
                    random_number2 = random.randint(-3, 3)
                    # Get the new size of the room
                    new_width = max(
                        min(
                            width + random_number1,
                            self.rooms_range[j][2],
                        ),
                        self.rooms_range[j][1],
                    )
                    new_depth = max(
                        min(
                            depth + random_number2,
                            self.rooms_range[j][4],
                        ),
                        self.rooms_range[j][3],
                    )
                    # Create a new room with the new size
                    new_room = (room_type, new_width, new_depth)
                    # replace the room with the new size
                    rooms[j] = new_room

            # 2. put it into Position_OnePlusOneEA to get the best layout for current sizes
            one_one_EA = Position_OnePlusOneEA(self.boundary, rooms)
            best_layout, fitness = one_one_EA.search()
            if fitness > self.best_fitness:
                self.best_layout = best_layout
                self.best_fitness = fitness
                self.best_layout_rooms = rooms.copy()
                # save snapshots every time better layout is found
                # force the count number to be like: 0000, 0025, 0050, ...
                save_snapshots(
                    current_time,
                    self.best_layout,
                    str(i).zfill(4),
                    self.best_fitness,
                    0,
                    0,
                    0,
                )
                print(
                    f"Iteration: {i}, Polygon Type: {self.best_layout.cluster.geom_type}, Fitness: {self.best_fitness}"
                )

            if i == self.max_iter - 1:
                save_snapshots(
                    current_time,
                    self.best_layout,
                    str(i).zfill(4),
                    self.best_fitness,
                    0,
                    0,
                    0,
                )
                print(
                    f"Iteration: {i}, Polygon Type: {self.best_layout.cluster.geom_type}, Fitness: {self.best_fitness}"
                )

        return self.best_layout, self.best_fitness


# List of configurations
# (mcts_iter, num_particles, pso_iter)
configurations = [
    # (1000, 200, 200),
    # (150, 200, 200),
    (150, 100, 2000),
    # (150, 2000, 2000),
    # (150, 600, 200),
    # (300, 20, 2000),
    # (300, 200, 200),
    # (300, 200, 2000),
    # (300, 400, 400),
]


def run_pso_mcts(mcts_iter, num_particles, pso_iter):
    print(
        f"Running PSO with 1+1EA_iter=100, num_particles={num_particles}, pso_iter={pso_iter}"
    )
    # algorithm logic
    # pso = PSO(
    #     rooms_range, num_particles=num_particles, max_iter=pso_iter, mcts_iter=mcts_iter
    # )
    # best_house = pso.optimize(boundary)

    # 1+1 EA
    size_one_one_EA = Size_OnePlusOneEA(boundary, rooms_range)
    best_house, fitness = size_one_one_EA.search()

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
