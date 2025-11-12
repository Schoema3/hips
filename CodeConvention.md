
The goal of this convention is to ensure consistent coding throughout the project.
This will improve readability and maintainability, as consistent styles and patterns reduce developers' cognitive load while they are working with the code.
Deviating from this code convention is acceptable, but only if absolutely necessary and unavoidable. 
# Coding guidelines: 
Apart from commentary and documentation, a developer can only express their intentions using three limiting tools.
These tools are the use of whitespace and indentation; the use of code keywords or programming language-specific concepts (e.g. the keyword const in C++ code); and the naming of variables.
Since written code is read by humans and machines alike, these three tools offer different degrees of freedom for expressing intentions.  

# White spaces and Indentation.
Unlike Python, C and C++ code compilers ignore white spaces and indentation, so there is a high degree of freedom to express oneself with these tools in code.
However, to achieve the main goal of these code conventions, which is to keep the code consistent, please adhere to the following rules. 
- Maximum number of characters per line: 80
- Indentation: [4 spaces] 
- Do not use trailing white spaces.
- & (reference) and * (pointer) are part of the type and not the name of the variable. Therefore, '*int& anInt*;' is okay, but '*int &anInt*;' is not.
- Don't overdue it with the indentation. 4 levels should be enough, every level more is a good indicator for some need to restructure. (https://youtu.be/CFRhGnuXG-4)
# Names for variables, objects and functions
Apart from some rules about special characters and syntax, naming things is where freedom of expression is at its highest. It is acceptable to use names consisting of a single character, such as the well-known for-loop *int i* to completeSentencesInOneSingeWordForAPublicFunction(). The main goal of naming things is to convey intention and meaning to the next person who reads it. As a rule of thumb, a name should be short if it is used often; if it is used rarely, the name should be long enough to make it clear what it does without having to read the documentation( or the code).

To use the syntax of names to transport more meaning, here are some rules to keep the code consist:

| What                        | Naming convention         | Example                        |
| --------------------------- | ------------------------- | ------------------------------ |
| Class names                 | PascalCase                | `EddyGenerator`                |
| Objects names               | lower_case_with_dashes    | `eddy_generator`               |
| Local variables             | lower_case_with_dashes    | `number_of_eddies`             |
| Members variables           | name starts with single _ | `_field_name`, `_checkField()` |
| Public Methods              | camelCase                 | `object.methodExample()`       |
| Non public Methods          | _camelCase                | `object._privateFunEx()`       |
| Boolean methods             | camelCase                 | `_isThis()` or `_hasThis()`    |
| Type indirection            | type* name or type& name  | `int* a` / `int& a`            |
| Template parameters         | Pascal_Case_With_Dashes   | `template< class Class_A >`    |
| Magic Numbers/Math constans | ALL_CAPS                  | `MIN_EDDY_SIZE`                |

Combinations of the above mentioned are possible. For example a class private magic number shall look like this: 
```
const int _MIN_EDDY_SIZE = 6; 
```

# Documentation and commentary inside the code 
The written documentation of classes, methods and functions are part of the header files. It shall only contain the description of  *what* the specific this is supposed to do ==not *how*== it is done. If the how needs explanation then this shall be done through a (block) commentary style in the source code. Sometimes it is helpful to reference papers with their DOI-link in the document, to indicate further reading material on what and how the specific thing should work.

The documentation inside the header files shall be compatible with doxygen like:
`///`
`/// ... text ...`
`///`
(ref.: https://www.doxygen.nl/manual/docblocks.html#cppblock)

Commentary inside the code is only helpful if it is used sparsely otherwise it is a good place to hide important information in a lot of clutter and noise, 
Therfore:
 - Don't leaf comments inside the code which explain language concepts or are trivial, like: 
   ```
    // Resize vectors to the number of variables
    varData.resize(nVar);
	```
	An exception to this would be regular expression, these always lack the ability of conveying meaning just by them self or when the trivial code needs a comment to convey the intention like this example form the python style guide [pep8](https://peps.python.org/pep-0008/#comments):
	```
	x = x + 1; // Compensate for the boarder
    ```
  

Magic numbers inside of code shall have commentary or documentation attached to them:
```
const int _MIN_EDDY_SIZE = 6; /// Minimal size of an Eddy 
```
