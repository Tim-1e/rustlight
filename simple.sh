cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s ualpha
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o ualpha.png -m 0.2 .\scene\meeting_ply.pbrt plane-single -n 40960 -s ualpha
cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o smis_jacobian_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s smis_jacobian -k 2 -x
cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s smis_all -k 2 -x
cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o cmis.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s cmis
cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o Average.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s average
cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s proxy_sample
cargo run --features="pbrt openexr" --release -- -t -2 -n 1 -o smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 40960 -s smis_all -k 2 -x
proxy Elapsed Integrator: 32667 ms
sims Elapsed Integrator: 18328 ms

plus
cargo run --features="pbrt openexr" --release -- -t -5 -n 2 -o ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 81920 -s ualpha
cargo run --features="pbrt openexr" --release -- -t -5 -n 10 -o ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane_single -n 409600 -s ualpha