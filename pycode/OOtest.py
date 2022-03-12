class vect:
    def __init__(self, x,y,z):
        self.x = x
        self.y = y
        self.z = z
       
    def ConstMult(self, const):
        temp=vect(0,0,0)
        temp.x=self.x*const
        temp.y=self.y*const
        temp.z=self.z*const
        return temp
    
    def dot(self, vec):
        temp = self.x * vec.x +self.y * vec.y + self.z * vec.z
        return temp
 
    def __repr__(self):
        return f" x-component: {self.x}, y-component: {self.y}, z-component: {self.z}"




v1 = vect(1,3,5)
v2 = vect(-1,-3,-5)
const = 10


print(v1)
print(v2)
print(v1.x, v1.y, v1.z)
temp=v1.ConstMult(const)
print(temp)
v1dotv2 =v1.dot(v2)
print (v1dotv2)