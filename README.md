# FastContainter
The cases when you need to get the nearest value of some array to the input variable are usual things. For this type of task the FastContainer was developed.
If you have a vector with values and its size is N, to solve this task it spents N actions per one searching O(N).
But what if the solution can be found almsot by one step? For example, std::unorder_map looks like a solution. Actually yes, but no.
Unordered map can evaluate the value only for the exact input: it doesn't snap the key value to the closest one.
FastContainer divide full range of values into several batches of fixed size to search out inside of this small batches. Finally, we have O(1*batch_size).
