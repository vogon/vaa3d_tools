#ifndef SEARCH_KDTREE_H
#define SEARCH_KDTREE_H

template <typename T>
class SearchKDTree
{
public:
	SearchKDTree(QList<T> & nodes, double (*location_fn)(const T &, int), 
		double (*d2_fn)(const T &, const T &)): location_fn(location_fn), d2_fn(d2_fn)
	{
		this->root = buildSubtree(nodes, 0, location_fn);
	}

	~SearchKDTree()
	{
		if (this->root != NULL) delete this->root;
	}

	T &nearestNeighbor(const T & pt) const
	{
		double r2 = std::numeric_limits<double>::infinity();
		return nnSearchSubtree(pt, this->root, r2, NULL)->value;
	}

	QList<T> rangeSearch(const T & center, double r2) const
	{
		QList<T> results;
		rangeSearchSubtree(center, this->root, r2, results);

		return results;
	}

private:
	class LessThanFunctor 
	{
	public:
		LessThanFunctor(double (*fn)(const T &, int), int axis): fn(fn), axis(axis) {}

		bool operator()(const T & left, const T & right) 
		{
			return this->fn(left, this->axis) < this->fn(right, this->axis);
		}

	private:
		double (*fn)(const T &, int);
		int axis;
	};

	struct Node 
	{
		double location;
		int axis;
		T value;
		Node *left, *right;

		Node() {}
		~Node() {
			if (this->left != NULL) delete this->left;
			if (this->right != NULL) delete this->right;
		}
	};

	Node *buildSubtree(QList<T> & nodes, int root_axis, double (*location_fn)(const T &, int))
	{
		// find median on the root axis
		QList<T> local_nodes(nodes);
		LessThanFunctor lt(location_fn, root_axis);

		qSort(local_nodes.begin(), local_nodes.end(), lt);

		int median_index = local_nodes.size() / 2;

		// scan left so that all the nodes equal to the median (if any) are to its right
		while (median_index > 0 && 
			(location_fn(local_nodes.at(median_index - 1), root_axis) ==
			 location_fn(local_nodes.at(median_index), root_axis)))
		{
			median_index--;
		}

		// create root node
		Node *root = new Node();
		root->location = location_fn(local_nodes.at(median_index), root_axis);
		root->axis = root_axis;
		root->value = local_nodes.at(median_index);

		// create child nodes
		int next_axis = (root_axis + 1) % 3;
		QList<T> left_nodes = local_nodes.mid(0, median_index);
		QList<T> right_nodes = local_nodes.mid(median_index + 1);

		if (!left_nodes.isEmpty())
		{
			root->left = buildSubtree(left_nodes, next_axis, location_fn);
		} 
		else 
		{
			root->left = NULL;
		}

		if (!right_nodes.isEmpty())
		{
			root->right = buildSubtree(right_nodes, next_axis, location_fn);
		} 
		else 
		{
			root->right = NULL;
		}

		return root;
	}

	void rangeSearchSubtree(const T & center, Node *root, double r2, QList<T> & results) const
	{
		if (root == NULL)
		{
			return;
		}

		// check to see if the root is closer to center than the cutoff
		double root_d2 = this->d2_fn(center, root->value);

		if (root_d2 < r2) {
			results.append(root->value);
		}

		// search the subtree of root that contains center
		if (this->location_fn(center, root->axis) < root->location)
		{
			// center is left of the root
			rangeSearchSubtree(center, root->left, r2, results);

			// if the r2-wide sphere intersects with the plane dividing the two subtrees,
			// search the other subtree as well
			double r = root->location - this->location_fn(center, root->axis);

			if ((r * r) < r2)
			{
				rangeSearchSubtree(center, root->right, r2, results);
			}
		}
		else
		{
			// pt is right of the root
			rangeSearchSubtree(center, root->right, r2, results);

			// if the r2-wide sphere intersects with the plane dividing the two subtrees,
			// search the other subtree as well
			double r = this->location_fn(center, root->axis) - root->location;

			if ((r * r) < r2)
			{
				rangeSearchSubtree(center, root->left, r2, results);
			}
		}
	}

	Node *nnSearchSubtree(const T & pt, Node *root, double & cutoff_r2, Node *default_closest) const
	{
		if (root == NULL)
		{
			return default_closest;
		}

		Node *closest = default_closest;

		// check to see if the root is closer to pt than the cutoff
		double root_d2 = this->d2_fn(pt, root->value);

		if (root_d2 < cutoff_r2)
		{
			cutoff_r2 = root_d2;
			closest = root;
		}

		// search the subtree of root that contains pt
		if (this->location_fn(pt, root->axis) < root->location)
		{
			// pt is left of the root
			closest = nnSearchSubtree(pt, root->left, cutoff_r2, closest);

			// if the cutoff_r2-wide sphere intersects with the plane dividing the two subtrees,
			// search the other subtree as well
			double r = root->location - this->location_fn(pt, root->axis);

			if ((r * r) <= cutoff_r2)
			{
				closest = nnSearchSubtree(pt, root->right, cutoff_r2, closest);
			}
		}
		else
		{
			// pt is right of the root
			closest = nnSearchSubtree(pt, root->right, cutoff_r2, closest);

			// if the cutoff_r2-wide sphere intersects with the plane dividing the two subtrees,
			// search the other subtree as well
			double r = this->location_fn(pt, root->axis) - root->location;

			if ((r * r) <= cutoff_r2)
			{
				closest = nnSearchSubtree(pt, root->left, cutoff_r2, closest);
			}
		}

		return closest;
	}

	Node *root;
	double (*location_fn)(const T &, int);
	double (*d2_fn)(const T &, const T &);
};

#endif /* !SEARCH_KDTREE_H */