use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use rand::Rng;
use rand::prelude::SliceRandom;
use rand::prelude::IteratorRandom; 

// Calculate Jaccard similarity
fn jaccard_similarity<T: Eq + Clone + std::hash::Hash>(
    set1: &HashSet<T>,
    set2: &HashSet<T>,
) -> f64 {
    //Defined as the size of the intersection divided by the size of the union 
    let intersection_size = set1.intersection(&set2).count() as f64;
    let union_size = set1.union(&set2).count() as f64;
    if union_size == 0.0 {
        0.0 
    } else {
        intersection_size / union_size
    }
}

#[allow(unused_assignments)]
fn clustering_coefficient(node: usize, graph: &HashMap<usize, HashSet<usize>>) -> f64 {
    let binding = HashSet::new();
    let neighbors = graph.get(&node).unwrap_or(&binding);
    let neighbors_num = neighbors.len();

    // If fewer than 2 neighbors, return undefined (0.0)
    if neighbors_num < 2 {
        return 0.0;
    }

    let mut triangles_num = 0;

    // Iterate over pairs of neighbors to find triangles
    for &a in neighbors {
        for &b in neighbors {
            if a < b && graph.get(&a).unwrap_or(&HashSet::new()).contains(&b) {
                triangles_num += 1;
            }
        }
    }

    // Count the total number of pairs among the neighbors
    let possible_pairs_num = (neighbors_num * (neighbors_num - 1)) / 2;

    // Calculate clustering coefficient
    if possible_pairs_num > 0 {
        triangles_num as f64 / possible_pairs_num as f64
    } else {
        0.0 // if there are no possible pairs, clustering coefficient is 0
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Open file
    let file_path = "/Users/adrianne/Downloads/facebook_combined.txt";
    let file = File::open(&file_path)?;
    let reader = BufReader::new(file);

    let mut graph: HashMap<usize, HashSet<usize>> = HashMap::new();

    // Parse the input file to construct the graph
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<usize> = line
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();

        let node1 = parts[0];
        let node2 = parts[1];

        // Add edges to the graph
        graph.entry(node1).or_insert_with(HashSet::new).insert(node2);
        graph.entry(node2).or_insert_with(HashSet::new).insert(node1);
    }

    // Calculate the proportion of friends' friends that are also your friends
    let mut total_friends_of_friends = 0;
    let mut total_mutual_friends_of_friends = 0;

    // Iterate through each friend of the user
    for (_, connections) in &graph {
        for &friend in connections {
            // Check if any of their friends are also friends with the user
            if let Some(friends) = graph.get(&friend) {
                for &fof in friends {
                    if connections.contains(&fof) {
                        total_mutual_friends_of_friends += 1;
                    }
                    total_friends_of_friends += 1;
                }
            }
        }
    }

    let proportion = total_mutual_friends_of_friends as f64 / total_friends_of_friends as f64;
    println!("Proportion of friends' friends that are also your friends: {}", proportion);

    // Find two vertices with the most similar or most dissimilar sets of connections
    let mut max_similarity = 0.0;
    let mut min_similarity = 1.0;
    let mut most_similar = (String::new(), String::new());
    let mut most_dissimilar = (String::new(), String::new());

    for (node1, connections1) in &graph {
        for (node2, connections2) in &graph {
            if node1 != node2 {
                let similarity = jaccard_similarity(&connections1, &connections2);
            if similarity > max_similarity {
                max_similarity = similarity;
                most_similar = (node1.to_string(), node2.to_string()); // Store owned strings
            }
            if similarity < min_similarity {
                min_similarity = similarity;
                most_dissimilar = (node1.to_string(), node2.to_string()); // Store owned strings
            }
        }
    }
}

    println!("Most similar: {:?} (Similarity: {})", most_similar, max_similarity);
    println!("Most dissimilar: {:?} (Similarity: {})", most_dissimilar, min_similarity);

    println!("Clustering Coefficient:");
    for (node, _) in &graph {
        let clustering_coef = clustering_coefficient(*node, &graph);
        println!("Node{}:{}", node, clustering_coef);
    }
    println!("Random Walk:");

    let mut current_vertex = *graph.keys().next().unwrap();
    let mut random_walk_sequence = Vec::new(); // Store the sequence of vertices visited

    for _ in 0..10 {
        println!("Current vertex: {}", current_vertex);
        random_walk_sequence.push(current_vertex);
        current_vertex = random_walk_step(&graph, current_vertex);
    
        // Calculate path length after each step
        let path_length = random_walk_sequence.len();
        println!("Path Length: {}", path_length);
    }

    // Node importance
    let mut visit_counts: HashMap<usize, usize> = HashMap::new();
    for &vertex in &random_walk_sequence {
        *visit_counts.entry(vertex).or_insert(0) += 1;
    }

        // Print node importance
    for (vertex, count) in visit_counts.iter() {
        println!("Node Importance: Vertex {}: Visits {}", vertex, count);
    }
    
    Ok(())    
}

// Random walk
pub fn random_walk_step(graph: &HashMap<usize, HashSet<usize>>, current_vertex: usize) -> usize {
    let mut rng = rand::thread_rng();
    let outgoing_edges = match graph.get(&current_vertex) {
        Some(edges) => edges,
        None => return current_vertex, // Return the same vertex if no outgoing edges
    };

    if outgoing_edges.is_empty() || rng.gen::<f64>() < 0.1 {
        let all_vertices: Vec<_> = graph.keys().collect();
        **all_vertices.choose(&mut rng).unwrap() // Dereference twice to get usize
    } else {
        if rng.gen_bool(0.9) {
            *outgoing_edges.iter().choose(&mut rng).unwrap() // Choose a random neighbor and return it
        } else {
            let all_vertices: Vec<_> = graph.keys().collect();
            **all_vertices.choose(&mut rng).unwrap() // Dereference twice to get usize
        }
    }
}

//test part
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_jaccard_similarity() {
        let set1: HashSet<usize> = [1, 2, 3].iter().cloned().collect();
        let set2: HashSet<usize> = [2, 3, 4].iter().cloned().collect();
        assert_eq!(jaccard_similarity(&set1, &set2), 0.5);
    }
    #[test]
    fn test_clustering_coefficient() {
        // Create a graph with triangular relationships
        let mut graph: HashMap<usize, HashSet<usize>> = HashMap::new();
        graph.insert(1, [2, 3].iter().cloned().collect());
        graph.insert(2, [1, 3].iter().cloned().collect());
        graph.insert(3, [1, 2].iter().cloned().collect());

        // Calculate clustering coefficient for each node
        let clustering_coef_1 = clustering_coefficient(1, &graph);
        let clustering_coef_2 = clustering_coefficient(2, &graph);
        let clustering_coef_3 = clustering_coefficient(3, &graph);
        
        // Expected clustering coefficient for nodes in a triangular graph is 1.0
        assert_eq!(clustering_coef_1, 1.0);
        assert_eq!(clustering_coef_2, 1.0);
        assert_eq!(clustering_coef_3, 1.0);
    }

    #[test]
    fn test_random_walk_step() {
        let mut graph: HashMap<usize, HashSet<usize>> = HashMap::new();
        graph.insert(1, [2, 3].iter().cloned().collect());
        graph.insert(2, [1, 3].iter().cloned().collect());
        graph.insert(3, [1, 2].iter().cloned().collect());

        let mut visited_vertices = HashSet::new();
        let mut current_vertex = 1;
        for _ in 0..10 {
            visited_vertices.insert(current_vertex);
            current_vertex = random_walk_step(&graph, current_vertex);
        }
        // Ensure that the random walk visited multiple vertices
        assert!(visited_vertices.len() > 1);
    }
    #[test]
    fn test_node_importance() {
        // Simulate a random walk sequence
        let random_walk_sequence = vec![1, 2, 3, 2, 1, 4, 5, 1, 2, 3, 6, 1, 7, 8, 2];
        
        // Calculate node importance
        let mut visit_counts: HashMap<usize, usize> = HashMap::new();
        for &vertex in &random_walk_sequence {
            *visit_counts.entry(vertex).or_insert(0) += 1;
        }
        
        // Verify the node importance counts
        assert_eq!(visit_counts[&1], 4);
        assert_eq!(visit_counts[&2], 4);
        // Add more assertions for other vertices as needed
    }
    #[test]
    fn test_path_length() {
        // Simulate a random walk sequence
        let random_walk_sequence = vec![1, 2, 3, 2, 1, 4, 5, 1, 2, 3, 6, 1, 7, 8, 2];

        // Verify the path length
        assert_eq!(random_walk_sequence.len(), 15);
    }
} 