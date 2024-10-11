from enum import Enum

class Sections(Enum):
    PUBLIC_SPACES = 1
    # EXAMPLES: Living Room, Dining Room, Kitchen

    PRIVATE_SPACES = 2
    # EXAMPLES: Bedroom, Master Bedroom, Children's Room, Toilet, Bathroom

    SERVICE_SPACES = 3
    # EXAMPLES: Laundry, Storage, Garage

    TRANSITIONAL_SPACES = 4
    # EXAMPLES: Hallway, Foyer, Terrace

    OTHERS = 5
    # EXAMPLES: Study, Game Room, Gym, Home Theater, Balcony/Patio, Courtyard