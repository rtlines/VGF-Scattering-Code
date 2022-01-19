
import json

with open("test_dicionary.json",'r') as fp: 
    variable=fp.read()
    classmates=json.loads(variable)


print(classmates['First'],classmates['Last'],classmates['Grades'])


age1=classmates['Age']


print(age1)
