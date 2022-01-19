import json

AList1=[1,2,3]
AList2=[4,5,6]

classmates={'First':'Christine', 'Last':'Whitney','Age':25,'Birthday':'March 10','Grades':AList1}
print(classmates['First'],classmates['Last'],classmates['Grades'])


with open("test_dicionary.json",'w') as fp: 
    variable=json.dumps(classmates, indent=4)
    fp.write(variable)
