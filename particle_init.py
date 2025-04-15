import json
import random

def generate_particles_absolute_center(
    rows,
    cols,
    spacing,                  # Absolute spacing in simulation units
    region_center,            # Absolute position of the region center
    velocity=(0.0, 0.0),
    mass=1.0,
    domain_min=(0.0, 0.0),
    domain_max=(1.0, 1.0),
    jitter=0.0,
    output_file="particles.json"
):
    particles = []

    region_w = spacing * (cols - 1)
    region_h = spacing * (rows - 1)

    center_x, center_y = region_center
    start_x = center_x - region_w / 2
    start_y = center_y - region_h / 2

    # Check: is the region inside the domain?
    if not (domain_min[0] <= start_x <= domain_max[0] and
            domain_min[1] <= start_y <= domain_max[1]):
        print("⚠️ WARNING: Region starts outside domain!")

    pid = 0
    for i in range(rows):
        for j in range(cols):
            x = start_x + j * spacing
            y = start_y + i * spacing

            if jitter > 0.0:
                x += random.uniform(-jitter, jitter)
                y += random.uniform(-jitter, jitter)

            particle = {
                "id": pid,
                "position": [x, y],
                "velocity": list(velocity),
                "mass": mass
            }
            particles.append(particle)
            pid += 1

    data = {
        "domain": {
            "min": list(domain_min),
            "max": list(domain_max)
        },
        "particles": particles
    }

    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"✅ Saved {len(particles)} particles to {output_file}")
    print(f"📏 Region center: {region_center}")
    print(f"📏 Region size: width = {region_w:.3f}, height = {region_h:.3f}")


# ✅ Example usage: particles in the upper middle of a 1.0 x 1.0 domain
generate_particles_absolute_center(
    rows=3,                     # Number of rows(particles)
    cols=4,                     # Number of columns(particles)
    spacing=0.05,               # Spacing between particles (in domain units)
    region_center=(0.5, 0.8),   # Center of particle region (in domain units)
    velocity=(0.0, 0.0),        # Initial velocity of particles (units per second in domain space)
    mass=1.0,                   # Mass of each particle
    domain_min=(0.0, 0.0),      # Bottom-left corner
    domain_max=(1.0, 1.0),      # Top-right corner
    jitter=0.0,               # Optional random offset added to particle positions
    output_file="particles.json"
)