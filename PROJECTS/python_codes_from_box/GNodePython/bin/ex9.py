class Organism(object):
    def __init__(self, species = 'unknown'):
        self.species = species
    def __str__(self):
        print 'I am an organism of species %s' % self.species

class Mammal(Organism):
    def __init__(self):
        self.warmblooded = True
    def __str__(self):
        print 'I am a mammal of species %s' % self.species
        
        
class Dog(Mammal):
    def __init__(self, sound = 'bark'):
        self.sound = sound
        self.species = 'Canis'
    def __str__(self):
        print 'I am a dog and I say %s' % self.sound
        
class Cat(Mammal):
    def __init__(self, sound = 'meow'):
        self.sound = sound
        self.species = 'Felis Minor'
    def __str__(self):
        print 'I am a cat and I say %s' % self.sound

if __name__ == "__main__":
    misty = Cat()
    print '%s' % misty.species
    print misty.sound
    blackie = Cat('bark')
    print blackie.sound
    print dir(misty)
    print blackie.__str__()
    print 'hello world'
    garfield = Organism('cat')
    print garfield.__str__()
    
