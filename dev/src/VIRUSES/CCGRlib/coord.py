class Coord:
    def __init__(self, x=0, y=0):
        self.X = x
        self.Y = y

    def get_coords(self):
        return tuple([self.X, self.Y])