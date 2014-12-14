#ifndef AST_CODE_HEADER_GUARD
#define AST_CODE_HEADER_GUARD
#include <vector>
//AST code
class node
{
public:
	virtual int calculate(std::vector<int>& edgeValues) = 0;
	virtual ~node()
	{}
};
class binaryNode : public node
{
public:
	binaryNode(node* left, node* right)
	:left(left), right(right)
	{}
	virtual ~binaryNode()
	{
		delete left;
		delete right;
	}
protected:
	node* left, *right;
};
class plusNode : public binaryNode
{
public:
	plusNode(node* left, node* right)
	: binaryNode(left, right)
	{}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return left->calculate(edgeValues) + right->calculate(edgeValues);
	}
};
class timesNode : public binaryNode
{
public:
	timesNode(node* left, node* right)
	:binaryNode(left, right)
	{}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return left->calculate(edgeValues) * right->calculate(edgeValues);
	}
};
class divideNode : public binaryNode
{
public:
	divideNode(node* left, node* right)
	:binaryNode(left, right)
	{}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return left->calculate(edgeValues) / right->calculate(edgeValues);
	}
};
class minusNode : public binaryNode
{
public:
	minusNode(node* left, node* right)
	: binaryNode(left, right)
	{}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return left->calculate(edgeValues) - right->calculate(edgeValues);
	}
};
class negateNode : public node
{
public:
	negateNode(node* negated)
	:negated(negated)
	{}
	virtual ~negateNode()
	{
		delete negated;
	}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return -negated->calculate(edgeValues);
	}
private:
	node* negated;
};
class integerNode : public node
{
public:
	integerNode(int value)
	:value(value)
	{}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return value;
	}
private:
	int value;
};
class edgeNode : public node
{
public:
	edgeNode(int edge)
	:edge(edge)
	{}
	virtual int calculate(std::vector<int>& edgeValues)
	{
		return edgeValues[edge];
	}
private:
	int edge;
};
#endif
