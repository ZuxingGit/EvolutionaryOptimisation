from enum import Enum

class RoomType(Enum):
    LIVING_ROOM = 1          # Living Room
    DINING_ROOM = 2          # Dining Room
    KITCHEN = 3              # Kitchen
    BEDROOM = 4               # Bedroom
    MASTER_BEDROOM = 5       # Master Bedroom
    CHILDRENS_ROOM = 6       # Children's Room
    BATHROOM = 7              # Bathroom
    MASTER_BATHROOM = 8      # Master Bathroom
    GUEST_BATHROOM = 9       # Guest Bathroom
    STUDY = 10                # Study/Office
    GAME_ROOM = 11           # Game Room
    GYM = 12                  # Gym
    HOME_THEATER = 13        # Home Theater
    LAUNDRY_ROOM = 14        # Laundry Room
    STORAGE_ROOM = 15        # Storage Room
    GARAGE = 16              # Garage
    FOYER = 17                # Foyer
    HALLWAY = 18             # Hallway
    TERRACE = 19             # Terrace/Balcony
    FAMILY_ROOM = 20         # Family Room
    GUEST_ROOM = 21          # Guest Room
    UTILITY_ROOM = 22        # Utility Room
    WALK_IN_CLOSET = 23      # Walk-in Closet
    COURTYARD = 24           # Courtyard (Outdoor or central area)
    TOILET = 25              # Toilet

    def __str__(self):
        return self.name.replace('_', ' ')

# Example usage
def describe_room(room_type):
    return f"This is a {room_type.name.replace('_', ' ').lower()}, number: {room_type.value}"

# Test
# for room in roomType:
#     print(describe_room(room))
